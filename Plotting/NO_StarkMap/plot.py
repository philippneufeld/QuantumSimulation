# Philipp Neufeld, 2021-2022

import os

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    filename = "NOStarkMap_20220407-114437_ludwigsburg.h5"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/"

    fig = plt.figure()
    ax = fig.add_subplot(111)

    with h5py.File(os.path.join(dir_path, filename)) as file:
        attrs = dict(file.attrs)
        keys = list(file.keys())

        print(attrs)
        for key in [k for k in keys if k.isdecimal()]:
            data = file[key]

            energies = np.array(data["Energies"][::])
            eField = np.ones_like(energies) * float(data.attrs["Electric_Field"])
            
            ax.plot(eField, energies, "C0.", ms=1.5)

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))

    plt.show()
