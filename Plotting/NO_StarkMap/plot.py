# Philipp Neufeld, 2021-2022

import os
import argparse
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt

def unitVec(i, n):
    vec = np.zeros(n)
    vec[i] = 1
    return vec

if __name__ == '__main__':

    filename = "NOStarkMap_20220407-162632_calcc"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/"

    fig = plt.figure()
    ax = fig.add_subplot(111)

    with h5py.File(os.path.join(dir_path, filename + ".h5")) as file:
        attrs = dict(file.attrs)
        keys = list(file.keys())

        basis = file["basis"][::]
        maxR = np.max([b[2] for b in basis])
        basisR = np.zeros((int(maxR)+1, len(basis)))
        for i in range(len(basis)):
            basisR[int(basis[i, 2]), i] = 1

        print(attrs)
        for key in tqdm.tqdm([k for k in keys if k.isdecimal()]):
            try:
                data = file[key]
                energies = np.array(data["Energies"][::])
                states = np.array(data["States"][::]).T
                eField = np.ones_like(energies) * float(data.attrs["Electric_Field"])

                colors = ["C0" for e in energies]
                for i, ev in enumerate(states):
                    tmp = basisR * ev[None, :]
                    col_id = np.argmax(np.diag(np.matmul(tmp, tmp.T)))
                    colors[i] = f"C{col_id}"

                ax.scatter(eField, energies, 0.5, colors)
            except Exception as e:
                print(e)

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))

    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

    plt.show()
