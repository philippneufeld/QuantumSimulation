# Philipp Neufeld, 2021-2022
import sys

import subprocess


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


if __name__ == "__main__":

    filename = "NORydExTD_20220907-233123_calcc"
    dir_path = (
        "/home/PI5/pneufeld/remote_home/Masterarbeit/07_TimeDependence/01_Ion_modG"
    )
    path = os.path.join(dir_path, filename + ".h5")

    trajs = sorted([40, 80, 120, 260])

    with h5py.File(path) as file:
        fig1, (ax1, ax2) = plt.subplots(
            1, 2,
            gridspec_kw={"width_ratios": [1, 2.5]},
        )

        n = len(file.keys())
        freqs = np.zeros(n)
        pops = np.zeros_like(freqs)
        
        for k, v in file.items():
            idx = int(k)
            freq = v.attrs["frequency"][0]
            freqs[idx] = freq
            tmin = np.ceil(v["t"][-1] / 2 * freq) / freq
            tmax = v["t"][-1]
            rng = np.logical_and(v["t"][:] >= tmin, v["t"][:] <= tmax)
            pops[idx] = np.max(v["populations"][:][rng])

        # freqs = file["frequencies"][:]
        # pops = file["populations"][:]
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
                label=f"{freq_pprint(float(data.attrs['frequency']))}",
            )

        ax2.set_xlabel("Evolution time t (ms)")
        ax2.set_ylabel("Relative ion population ($10^{-7}$)")
        ax2.legend(ncol=5, loc="upper center")

    plt.show()
    
    # fig1.align_xlabels()
    # fig1.tight_layout()
    # fig1.savefig("fig1.pdf")
    # subprocess.run(["convert", "-density", "300", "fig1.pdf", "fig1.png"])
