# Philipp Neufeld, 2021-2022

import numpy as np
import scipy.constants as cs
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from tqdm import tqdm
from qsim.ladder import LadderSystem

def calc_doppler_ss(det, sys):
    ss = sys.steadystate([0.0, 0.0, det, 0.0])
    return np.imag(ss[0, 1])

if __name__ == '__main__':

    c = cs.speed_of_light

    lvlX = 0
    lvlA = lvlX + c / 226.97e-9
    lvlH = lvlA + c / 540e-9
    lvlR = lvlH + c / 834.92e-9
    lvlI = 9.27 * cs.e / cs.h

    rabiXA = 3932230.4421452461
    rabiAH = 11025372.503863104
    rabiHR = 55126862.519315511

    decayAX = 13.8e6
    decayHA = 10.0e6
    decayRH = 100.0e6
    decayRI = 6e5
    decayT = 407744.26954339806

    levels = [lvlX, lvlA, lvlH, lvlR, lvlI]
    rabis = [(rabiXA, 1.0), (rabiAH, -1.0), (rabiHR, -1.0), (0.0, -1.0)]
    decays = [
        (1, 0, decayAX + decayT),
        (2, 1, decayHA),
        (3, 2, decayRH),
        (3, 4, decayRI),
        (2, 0, decayT),
        (3, 0, decayT),
        (4, 0, decayT),
    ]
    mass = 4.9826301286306257e-26
    system = LadderSystem(levels, rabis, decays, mass)

    ex = ProcessPoolExecutor()
    dets = np.linspace(-500e6, 500e6, 501)
    pops = np.array([*tqdm(ex.map(calc_doppler_ss, dets, repeat(system)), total=len(dets))])

    plt.plot(dets, pops)
    plt.show()
