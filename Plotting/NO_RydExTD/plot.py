# Philipp Neufeld, 2021-2022
# Fabian Munkes, 2022 (did the most work) as always!
import sys

sys.path.append(
    "/mnt/ceph/file/groups/MicCells/TraceGasSensing/Measurements/2021/2021-Electric-Map-Eval"
)

import niceplot

niceplot.thesis_defaults()
niceplot.use_pgf()

import subprocess


from typing import Any, Tuple

import os
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt


def freq_pprint(val):
    units = {1: "Hz", 1e3: "kHz", 1e6: "MHz", 1e9: "GHz"}
    units = {
        1: r"\hertz",
        1e3: r"\kilo\hertz",
        1e6: r"\mega\hertz",
        1e9: r"\giga\hertz",
    }
    for v in reversed(sorted(units)):
        if val > v:
            # return f"${val / v:.2f}$ {units[v]}"
            return f"\\SI{{{val/v:.2f}}}{{{units[v]}}}"
    return "ERROR"


if __name__ == "__main__":

    # filename = "Sim_Fine_03"
    filename = "NORydExTD_20220503-093615_calcc"
    dir_path = (
        "/home/PI5/pneufeld/remote_home/Masterarbeit/07_TimeDependence/01_Ion_modG"
    )
    dir_path = (
        "/mnt/ceph/file/users/pneufeld/Masterarbeit/07_TimeDependence/01_Ion_modG"
    )
    path = os.path.join(dir_path, filename + ".h5")

    # trajs = [0, 100, 200, 300, 400]
    trajs = sorted([40, 80, 120, 160])

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

        # fig1 = plt.figure(figsize=niceplot.set_size(600, ratio=1/1.61803398875/2))
        fig1, (ax1, ax2) = plt.subplots(
            1,
            2,
            gridspec_kw={"width_ratios": [1, 2.5]},
            figsize=niceplot.set_size(800, ratio=1 / 1.61803398875 / 2),
        )
        # ax2 = fig1.add_subplot(122) #, sharey = ax1)
        ax1.loglog()
        ax1.plot(freqs, np.array(pops) * 1e7)
        ax1.set_xlabel("Chop frequency of intermediate laser (Hz)")
        ax1.set_ylabel(r"Mean relative ion population ($10^{-7}$)")

        print("Plotting...")
        for data in tqdm.tqdm(
            [file[k] for k in file.keys() if k.isdecimal() and int(k) in trajs]
        ):
            plt.plot(
                data["t"][:] * 1e6,
                data["populations"][:] * 1e7,
                label=f"${freq_pprint(float(data.attrs['frequency']))}$",
            )

        ax2.set_xlabel(r"Evolution time t (\si{\micro\s})")
        ax2.set_ylabel("Relative ion population ($10^{-7}$)")
        ax2.set_xlim((0, 1e-5 * 1e6))
        ax2.set_ylim((-.4, 4.9))
        ax2.legend(ncol=5, loc="upper center")

    # plt.show()
    fig1.align_xlabels()
    fig1.tight_layout()
    fig1.savefig("fig1.pdf")

    subprocess.run(["convert", "-density", "300", "fig1.pdf", "fig1.png"])
