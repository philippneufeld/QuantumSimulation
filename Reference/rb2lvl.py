# Philipp Neufeld, 2021-2022

import numpy as np
import scipy.constants
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from tqdm import tqdm
from qsim.ladder import LadderSystem

def calc_doppler_ss(det, sys):
    ss = sys.doppler_steadystate((det,))
    return np.imag(ss[0, 1])

if __name__ == '__main__':

    levels = [0, scipy.constants.speed_of_light / 780.241e-9]
    rabis = [(3.5e6, 1.0),]
    decays = [(1, 0, 6.065e6),]
    mass = 1.44316060e-25
    system = LadderSystem(levels, rabis, decays, mass)

    ex = ProcessPoolExecutor()
    dets = np.linspace(-1e9, 1e9, 501)
    pops = np.array([*tqdm(ex.map(calc_doppler_ss, dets, repeat(system)), total=len(dets))])

    plt.plot(dets, pops)
    plt.show()
