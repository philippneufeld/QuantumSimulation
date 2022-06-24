# Philipp Neufeld, 2021-2022

from typing import Any, Tuple

import os
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from glob import glob

if __name__ == '__main__':

    filename = "NOStarkMap_20220428-140003_calcc"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/"

    # filename = "NOStarkMap_20220619-143956_monaco"
    # dir_path = "/home/pneufeld/git/QuantumSimulation/build/apps/NO_StarkMap"
    
    path = os.path.join(dir_path, filename + ".h5")

    paths = sorted(glob(f"{dir_path}/NOStarkMap*.h5"))[-4:]

    for path in paths:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title(path.split("/")[-1])

        # color_palette = np.array([col.to_rgb(f"C{i}") for i in range(10)])
        color_palette = np.array([col.to_rgb(c) for c in ("black", "blue", "red", "green", "purple", "orange", "grey")])
                
        with h5py.File(path) as file:
            for datagroup in tqdm.tqdm([file[k] for k in file.keys() if k.isdecimal()]):

                energies = np.array(datagroup["Energies"][::])
                eField = np.ones_like(energies) * float(datagroup.attrs["Electric_Field"])
                character = np.array(datagroup["Character"][::])

                # colors = color_palette[np.mod(character[:, 2].astype(int), len(color_palette)), :]
                colors = color_palette[np.clip(character[:, 2].astype(int), 0, len(color_palette)-1), :]
                
                # ax.scatter(eField, energies, 0.5, colors)
                ax.scatter(eField, energies, 0.75, colors)
        
        ax.set_xlim((0, 25))
        ax.set_ylim((-66.5, -61.5))
        ax.set_xlabel("Electric field (V / cm)")
        ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")       

        # stark_map = NOStarkMap(os.path.join(dir_path, filename + ".h5"))
        # stark_map.plot()

        fig.tight_layout()

        # fig.savefig(os.path.join(dir_path, filename + ".pdf"))
        # fig.savefig(os.path.join(dir_path, filename + ".png"))

    plt.show()
