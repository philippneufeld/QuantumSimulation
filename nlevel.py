""" This module provides several functions for simulating an N-Level-System using QUTIP """
# pylint: disable=invalid-name, unused-argument
from datetime import datetime
from functools import partial
from itertools import combinations
import multiprocessing as mp
from typing import Callable, Iterable, Tuple, Type

import numpy as np
from qutip import basis, expect, steadystate, qeye
import scipy.constants as const
from tqdm.contrib.concurrent import process_map


def H(detunings: dict, rabis: dict, pop_op: dict, trans_op: dict) -> None:
    """
    Hamiltonian. This function will be defined dynamically in 'define_hamiltonian'

    args:
        detunings: A dictionary holding the detunings per transition
        rabis: A dictionary holding the rabi frequencies per transition
        pop_op: The population operators of all states
        trans_op: The transition operators of all states
    """


def L(decays: dict, dec_op: dict) -> None:
    """Lindblad. This function will be defined dynamically in 'define_lindblad'

    args:
        decays: The decay values per transition (top to bottom)
        dec_op: The decay operators per transition (top to bottom)
    """


def get_states(names: Iterable[str]) -> dict:
    """
    Generates states kets of str in Iterable

    args:
        names: States description
    """

    num_states = len(names)
    states = {}
    for cnt, name in enumerate(names):
        states[name] = basis(num_states, cnt)
    return states


def get_population_operators(states: dict) -> dict:
    """
    Generates a dictionary of all population operators

    args:
        states: states dictionary
    """
    opop = {}
    for state, value in states.items():
        _op = value * value.dag()
        opop[state] = _op
    return opop


def get_transition_operators(states: dict, transitions: Iterable[Iterable]) -> dict:
    """
    Generates a dictionary of all transition operators

    args:
        states: states dictionary
        transitions: tuple of parameter class with (init state (str), final state (str), strength (float), k (float))
    """
    names = [f"{x}{y}" for x, y, _, _ in transitions]
    otrans = {}
    for name, (_si, _sf, _, _) in zip(names, transitions):
        si = states[_si]
        sf = states[_sf]
        t = (si * sf.dag() + sf * si.dag()) * 0.5
        otrans[name] = t

    return otrans


def get_transition_values(transitions: Iterable[Iterable]) -> dict:
    """Generates a dictionary where d[transition] = value

    args:
        transitions: see parameter class
    """
    names = [f"{x}{y}" for x, y, _, _ in transitions]
    t = {}
    for name, (_, _, v, _) in zip(names, transitions):
        t[name] = v
    return t


def get_decay_operators(states: dict, decays: Iterable[Iterable]) -> dict:
    """
    Generates a dictionary of all decay operators

    args:
        states: states dictionary
        decay: list of tuples of parameter class with (init state (str), final state (str), strength (float))
    """
    names = [f"{x}{y}" for x, y, _ in decays]
    odec = {}
    for name, (_si, _sf, _) in zip(names, decays):
        si = states[_si]
        sf = states[_sf]
        t = sf * si.dag()
        odec[name] = t

    return odec


def get_decay_values(decays: Iterable[Iterable]) -> dict:
    """Generates a dictionary where d[decay] = value

    args:
        decays: see parameter class
    """
    names = [f"{x}{y}" for x, y, _ in decays]
    d = {}
    for name, (_, _, v) in zip(names, decays):
        d[name] = v
    return d


def get_wavevectors(transitions: Iterable[Iterable]) -> dict:
    """Generates a dictionary of wavevectors where d[transition] = value

    args:
        transitions: see parameter class
    """
    names = [f"{x}{y}" for x, y, _, _ in transitions]
    k = {}
    for name, (_, _, _, v) in zip(names, transitions):
        k[name] = v
    return k


def define_hamiltonian(pop_op: dict, trans_op: dict) -> None:
    """This functions defines the Hamiltonian

    The Hamiltonian is defined dynamically based on all population and
    transition operators. Thus exec is used to parse the string. Since
    multiprocessing can only work with top level defined objects, the
    globals dictionary is expanded. See pickling.

    See the defintion of :func:`H` for a list of arguments to be passed to the generated
    Hamiltonian.

    args:
        pop_op: Dictionary of population operators
        trans_op: Dictionary of transition operators
    """

    s = "def H(detunings: dict, rabis: dict, pop_op: dict, trans_op: dict): return 2 * np.pi * ("
    # there are less transitions than operators
    for transition, population in zip(trans_op, pop_op):
        s += f" - detunings['{transition}'] * pop_op['{population}'] + rabis['{transition}'] * trans_op['{transition}']"

    s += ")"
    exec(s)
    globals()["H"] = locals()["H"]


def define_lindblad(decay_op):
    """This functions defines the Lindblad.

    See define_hamiltonian for explanation.

    See the defintion of :func:`L` for a list of arguments to be passed to the generated
    Hamiltonian.

    args:
        decay_op: Dictionary of decay operators
    """
    s = "def L(decays: dict, dec_op: dict): return ["
    for dec in decay_op:
        s += f"np.sqrt(decays['{dec}'] * 2 * np.pi) * dec_op['{dec}'],"
    s += "]"
    exec(s)
    globals()["L"] = locals()["L"]


def get_contributing_velocities(
    wavevectors: dict,
    detuning: float,
    velocities_full: np.ndarray,
    thermal_full: np.ndarray,
    decays: dict,
    width_factor: float = 100,
    check: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """This functions calculates the velocity intervals contributing to a
    specific detuning. A width dependend on the largest decay is used as
    range around the calculated velocity.

    args:
        wavevectors: Dictionary of all wavevectors
        detuning: Laser detuning of the scanned transition.
        velocities_full: The numpy array of the full velocity range simulated.
        thermal_full: The thermal function matching the full velocity range.
        decays: Dictionary of decays
        width_factor: Times the width around the calculated velocity
        check: Print's some nice stuff for debugging if turned on.
    returns:
        Tuple[np.ndarray, np.ndarray]: velocities and thermal range shrinked to relevant contributions
    """
    k = wavevectors
    num_wavevectors = len(wavevectors)
    # try to get contributing velocities:
    if num_wavevectors == 1:
        if check:
            print(f"Only one transition, contributing velocities as is")
        return [velocities_full], [thermal_full], ["blubb"]

    fields = list(wavevectors.keys())
    comb = []
    comb.append([fields[0]])
    for i in range(1, len(fields[1:]) + 1):
        for c in combinations(fields[1:], i):
            comb.append([fields[0], *c])
    # for i in range(1, len(fields) + 1):
    #     for c in combinations(fields, i):
    #         comb.append(c)
    if check:
        print("Contributing velocities:")
    v_dets = []
    v_str_contribs = []
    for c in comb:
        _sum_k = 0
        _s = []
        for _k in c:
            _sum_k += wavevectors[_k]
            _s.append(f"k_{_k}")
        v_str_contribs.append('+'.join(_s))
        v_dets.append(detuning / _sum_k)
        if check:
            print(f"detuning/({' + '.join(_s)})")

    b = width_factor
    c = decays["AX"] / k["XA"] * b
    _d = 0
    _k = 0
    for c in fields:
        _d += decays[c[::-1]]
        _k += k[c]
    c = _d / _k * b

    velocities, thermals = [], []
    for v_det in v_dets:
        _range = np.nonzero((velocities_full > v_det - c) & (velocities_full < v_det + c))
        velocities.append(velocities_full[_range])
        thermals.append(thermal_full[_range])

    if check:
        print(f"Contributing velocity steps: {velocities[0].size}")

    return velocities, thermals, v_str_contribs


def get_mean_velocity(mass: float, temperature: float) -> float:
    """
    Calculates the mean velocity of a MB distribution.

    args:
        mass: Mass in kg
        temperature: Temperature in °C
    returns:
        float: mean velocity
    """
    return np.sqrt(8 * const.k * (temperature + 273.15) / (mass * np.pi))


def get_mean_velocity_NO(temperature: float) -> float:
    """
    args:
        temperature: Temperature in °C
    returns:
        float: mean velocity of NO
    """

    return get_mean_velocity(4.9816e-26, temperature)


def get_thermal_dist_NO() -> Callable[[float, float, float], np.ndarray]:
    """Returns a function for calculating the thermal distribution of
    velocities for NO.

    returns:
        Callable[[float, float, float], np.ndarray]: Thermal function
    """

    def thermal(vel: np.ndarray, temperature: float, dv: float) -> np.ndarray:
        """
        args:
            vel: Array of velocities
            temperature: Temperature in °C
            dv: normalization factor
        returns:
            np.ndarray: thermal distribution
        """
        return (
            dv
            * np.exp(
                -(vel ** 2) * (0.001804 / (273.15 + temperature))
            )  # -mv^2 / (2*k_B T)
            * 0.023964  # sqrt(m/(k_B * 2pi))
            * np.sqrt(1 / (273.15 + temperature))
        )

    return thermal


class __Physics:
    """
    Helper class containing all properties for simulating, since everything has
    to be defined a toplevel for multiprocessing to work.
    """

    def __init__(self, parameters):
        """The init function allows to inherit class attributes (deepcopy does
        not help here). Thus we can copy the parameter class
        args:
            parameters: Parameter class
        """
        for p, v in parameters.__dict__.items():
            if not p.startswith("__"):
                self.__dict__[p] = v


def get_physics(PARAMETERS) -> Type[__Physics]:
    """Calculates all operators, decays etc. and fills helper class

    args:
        PARAMETERS: parameter class
    returns:
        Type[__Physics]: Parameter class expanded by operators and dictionaries of type d[transition] = value etc.
    """
    p = PARAMETERS
    physics = __Physics(p)
    physics.decays = get_decay_values(p.decays)
    physics.transitions = get_transition_values(p.transitions)

    physics.states = get_states(p.letters)

    physics.opop = get_population_operators(physics.states)
    physics.otrans = get_transition_operators(physics.states, p.transitions)
    physics.odec = get_decay_operators(physics.states, p.decays)

    define_hamiltonian(physics.opop, physics.otrans)
    define_lindblad(physics.odec)

    physics.k = get_wavevectors(p.transitions)
    physics.wavevectors = physics.k

    physics.mean_velocity = get_mean_velocity_NO(p.temperature)

    return physics


def calc(detuning, phys):
    scan_det = {}
    for transition in phys.transitions:
        scan_det[transition] = 0

    scan_det[f"{''.join(phys.scan_transition)}"] = detuning

    velocities, thermals, velo_contribs_desc = get_contributing_velocities(
        phys.k,
        detuning,
        phys.velocities,
        phys.thermal_dist,
        phys.decays,
        width_factor=phys.velocity_width_factor,
    )

    chi = 0
    populations = {}
    for state in phys.states:
        populations[state] = 0.0

    populations_per_contrib = {}
    for desc in velo_contribs_desc:
        populations_per_contrib[desc] = {}
        for state in phys.states:
            populations_per_contrib[desc][state] = 0.0

    for subvelos, subthermals, vdesc in zip(velocities, thermals, velo_contribs_desc):
        for speed, therm in zip(subvelos, subthermals):
            _detunings = {}
            for transition in phys.transitions:
                _detunings[transition] = speed * phys.k[transition] + scan_det[transition]
            curH = H(
                detunings=_detunings,
                rabis=phys.transitions,
                pop_op=phys.opop,
                trans_op=phys.otrans,
            )

            curL = L(decays=phys.decays, dec_op=phys.odec)

            ss = steadystate(
                curH,
                curL,
                method="direct",
            )

            n_states = len(phys.states)
            chi += (
                expect(
                    1 / n_states * qeye(n_states),
                    ss,
                )
                * therm
            )

            ### TODO: 2D plot
            exps = []
            for pop in populations:
                exp = expect(phys.opop[pop], ss) * therm
                exps.append(exp)
                populations[pop] += exp
            for pop, exp in zip(populations, exps):
                populations_per_contrib[vdesc][pop] += exp

    return {
        "chi": chi,
        "populations": populations,
        "populations_per_contrib": populations_per_contrib,
        "detuning": detuning,
        "velocitiesSize": velocities[0].size,
    }


def run_simulation(params):
    phys = get_physics(params)
    # debug
    get_contributing_velocities(
        phys.k,
        0.5,
        phys.velocities,
        phys.thermal_dist,
        phys.decays,
        width_factor=phys.velocity_width_factor,
        check=True,
    )

    start = datetime.now()
    results = process_map(
        partial(calc, phys=phys),
        phys.laser_detunings,
        max_workers=mp.cpu_count() - 1,
        chunksize=1,
    )

    lres = len(results)

    # sort by state, not iteration
    populations_final = {}
    populations_per_contrib_final = {}
    for desc in results[0]["populations_per_contrib"]:
        populations_per_contrib_final[desc] = {}
        for population in results[0]["populations"]:
            populations_per_contrib_final[desc][population] = np.empty(lres)

    for population in results[0]["populations"]:
        populations_final[population] = np.empty(lres)

    chi = np.empty(lres, dtype=np.cfloat)

    for c, result in enumerate(results):
        pop = result["populations_per_contrib"]
        for desc in pop:
            for state in populations_final:
                populations_per_contrib_final[desc][state][c] = pop[desc][state]

    for c, result in enumerate(results):
        pop = result["populations"]
        for state in populations_final:
            populations_final[state][c] = pop[state]
        chi[c] = result["chi"]
        velocitiesSize = result["velocitiesSize"]

    final_result = {
        "population": populations_final,
        "populations_per_contrib": populations_per_contrib_final,
        "chi": chi,
        "detunings": phys.laser_detunings,
        "velocitiesSize": velocitiesSize,
    }

    end = datetime.now()
    delta = end - start
    seconds = delta.seconds
    minutes = int(seconds / 60)
    seconds = seconds % 60
    print(f"Done in {minutes}min and {seconds}s")

    return final_result, end
