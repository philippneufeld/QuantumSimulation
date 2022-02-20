# Philipp Neufeld, 2021-2022

from typing import List, Tuple
import numpy as np
import scipy.constants
import scipy.integrate
import qutip



class LadderSystem:

    def __init__(self, levels: List[float], couplings: List[Tuple[int, int, float]], 
                 decays: List[Tuple[int, int, float]], mass: float, temperature: float=300) -> None:
        """Quantum ladder system
        Quantum system that has levels in a ladder-like configuration

        Args:
            dims (int): Number of levels in the syste
            couplings (List[Tuple[float, float]]): rabi frequencies and propagation direction of the couplings between the levels
            decays (List[Tuple[int, int, float]]): Decay rates. These decays shall be
                passed via a list that contains the start and end state of the decay as 
                well as the decay rate
            mass (float): Mass of the atom/molecule in kg
            temperature (float): Temperature of the system in K

        Raises:
            ValueError: Raised in case of an invalid (or inconsistent) argument
        """
        self.dims = len(levels)
        if any(levels[i+1] <= levels[i] for i in range(self.dims-1)):
            raise ValueError("Levels must be storted and non-degenerate")
        self.levels = np.array(levels)
        
        if len(couplings) != (self.dims - 1):
            raise ValueError("Invalid number of rabi frequencies")
        self.couplings = couplings

        if any(d[0] >= self.dims or d[1] >= self.dims or d[0]==d[1] for d in decays):
            raise ValueError("Invalid decays")
        self.decays = decays

        if (mass <= 0.0):
            raise ValueError("Mass must be positive")
        self.mass = mass

        if (temperature < 0.0):
            raise ValueError("Temperature may not be negative")
        self.temperature = temperature

    def operator_basis(self, i: int, j: int) -> qutip.Qobj:
        """Create a basis "vector" in the operator space

        Args:
            i (int): row of the nonzero element
            j (int): column of the nonzero element

        Returns:
            qutip.Qobj: calculated operator
        """
        return qutip.basis(self.dims, i) * qutip.basis(self.dims, j).dag()

    def hamiltonian(self, detunings: List[float]) -> qutip.Qobj:
        """Calculate the hamiltonian (in units of frequencies) for the given detunings

        Args:
            detunings (List[float]): Coupling laser detunings

        Raises:
            ValueError: raised if wrong number of detunings is passed

        Returns:
            qutip.Qobj: calculated hamiltonian
        """

        if len(detunings) != (self.dims - 1):
            raise ValueError("Invalid number of detunings")

        H = qutip.qzero(self.dims)

        # atom and light hamiltonian
        for i, el in enumerate(np.cumsum(detunings)):
            H -= el * self.operator_basis(i+1, i+1)

        # atom-light interaction hamiltonian
        for i, coupling in enumerate(self.couplings):
            rabi, _ = coupling
            H += 0.5 * rabi * self.operator_basis(i, i+1)
            H += 0.5 * np.conjugate(rabi) * self.operator_basis(i+1, i)

        return H * (2 * np.pi)

    def lindblad_operators(self) -> List[qutip.Qobj]:
        """Calculates the lindblad jump operators

        Returns:
            List[qutip.Qobj]: Linblad jump operators
        """
        ops = []
        for i, f, rate in self.decays:
            ops.append(np.sqrt(2*np.pi*rate) * self.operator_basis(f, i))
        return ops

    def steadystate(self, detunings: List[float]) -> qutip.Qobj:
        """Calculate the steady state of the system

        Args:
            detunings (List[float]): Coupling laser detunings

        Returns:
            qutip.Qobj: Steady state density matrix
        """
        H = self.hamiltonian(detunings)
        Ls = self.lindblad_operators()
        return qutip.steadystate(H, Ls, method="direct")        

    def detuning_doppler_shift(self, detunings: List[float], velocity: float) -> List[float]:
        """Calculate the doppler shifted detunings

        Args:
            detunings (List[float]): detunings to be doppler adjusted
            velocity (float): velocity of the frame of reference

        Raises:
            ValueError: Raised for invalid arguments

        Returns:
            List[float]: Adjusted detunings
        """
        
        if len(detunings) != (self.dims - 1):
            raise ValueError("Invalid number of detunings")

        # freq' = freq - k*v = freq * (1 +- v/c)
        # freq = splitting + detuning
        
        splittings = self.levels[1:] - self.levels[:-1]        
        freqs = np.array(detunings) + splittings
        factors = np.array([1 - coup[1]*velocity/scipy.constants.speed_of_light for coup in self.couplings])
        adjFreqs = freqs * factors

        return adjFreqs - splittings

    def doppler_steadystate(self, detunings: List[float], steps: int=1001) -> qutip.Qobj:
        """Calculate the doppler integrated steady state of the system

        Args:
            detunings (List[float]): Coupling laser detunings
            steps (int): Integration steps

        Returns:
            qutip.Qobj: Steady state density matrix
        """

        # should be an odd number
        steps = steps if steps % 2 == 1 else steps + 1
        
        nsigmas = 3.5
        sigma = np.sqrt(scipy.constants.Boltzmann*self.temperature/self.mass)
        vs = np.linspace(-nsigmas*sigma, nsigmas*sigma, steps)
        ws = 1 / (np.sqrt(2*np.pi)*sigma) * np.exp(-vs**2/(2*sigma**2))

        y = [ws[i]*self.steadystate(self.detuning_doppler_shift(detunings, vs[i])) for i in range(steps)]
        return np.sum(y, axis=0) * (2*nsigmas*sigma / (steps-1))



