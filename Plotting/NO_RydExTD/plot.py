# Philipp Neufeld, 2021-2022

from typing import Any, Tuple

import os
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    # filename = "Sim_Fine_03"
    filename = "NORydExTD_20220503-091745_ludwigsburg"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/07_TimeDependence/01_Ion_modG"
    path = os.path.join(dir_path, filename + ".h5")

    # trajs = [0, 100, 200, 300, 400]
    trajs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    with h5py.File(path) as file:
        
        freqs = file["frequencies"][:]
        pops = file["populations"][:]

        plt.figure()
        plt.loglog()
        plt.plot(freqs, pops)
        plt.xlabel("Chop frequency (green laser)")
        plt.ylabel("Mean relative ion population")
        
        plt.figure()
        for data in tqdm.tqdm([file[k] for k in file.keys() if k.isdecimal() and int(k) in trajs]):
            plt.plot(data["t"][:], data["populations"][:], label = f"$f = {float(data.attrs['frequency']):.3e}$")

        plt.xlabel("t")
        plt.ylabel("Relative ion population")
        plt.xlim((0, 2e-6))
        plt.legend()

    plt.show()
            

