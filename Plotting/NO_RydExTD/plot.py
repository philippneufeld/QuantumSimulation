# Philipp Neufeld, 2021-2022

from typing import Any, Tuple

import os
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt

def freq_pprint(val):
    units = {1: "Hz", 1e3: "kHz", 1e6: "MHz", 1e9: "GHz"}
    for v in reversed(sorted(units)):
        if val > v:
            return f"${val / v:.2f}$ {units[v]}"
    return "ERROR"

if __name__ == '__main__':

    # filename = "Sim_Fine_03"
    filename = "NORydExTD_20220503-093615_calcc"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/07_TimeDependence/01_Ion_modG"
    path = os.path.join(dir_path, filename + ".h5")

    # trajs = [0, 100, 200, 300, 400]
    trajs = [0, 40, 80, 120, 160]

    with h5py.File(path) as file:
        
        freqs = file["frequencies"][:]
        pops = file["populations"][:]

        # do averaging of trajectories
        # pops = np.zeros_like(freqs)
        # for i, data in tqdm.tqdm(enumerate([file[k] for k in sorted(file.keys()) if k.isdecimal()])):
        #     traj_p = np.array(data["populations"][:])
        #     traj_t = np.array(data["t"][:])
        #     tstart = 0.0 / freqs[i]
        #     rng = traj_t > tstart
        #     tsim = np.max(traj_t[rng]) - np.min(traj_t[rng])
        #     dt = traj_t[1]-traj_t[0]     
        #     pops[i] = traj_p[rng].sum() / tsim * dt

        plt.figure()
        plt.loglog()
        plt.plot(freqs, pops)
        plt.xlabel("Chop frequency (green laser)")
        plt.ylabel("Mean relative ion population")
        
        plt.figure()

        print("Plotting...")
        for data in tqdm.tqdm([file[k] for k in file.keys() if k.isdecimal() and int(k) in trajs]):
            plt.plot(data["t"][:], data["populations"][:], label = f"$f = {freq_pprint(float(data.attrs['frequency']))}$")

        plt.xlabel("t")
        plt.ylabel("Relative ion population")
        plt.xlim((0, 1e-5))
        plt.legend()

    plt.show()
            

