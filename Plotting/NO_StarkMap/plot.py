# Philipp Neufeld, 2021-2022

import os
import socket
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from glob import glob

COLOR_PALETTE = np.array([col.to_rgb(c) for c in ("black", "blue", "red", "green", "purple", "orange", "grey")])

def plot_starkmap_scatter(path):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efield = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            character = np.array(datagroup["Character"])

            colors = COLOR_PALETTE[np.clip(character[:, 2].astype(int), 0, len(COLOR_PALETTE)-1), :]
            ax.scatter(efield*np.ones(energies.shape), energies, 0.75, colors)

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

    return fig, ax

def plot_starkmap_lines(path, color_mode=2):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        
        efields = np.empty(len(keys))
        stark_map = []

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efields[idx] = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            characters = np.array(datagroup["Character"])

            if idx == 0:
                stark_map = [[(c, [e])] for c, e in zip(characters[:, color_mode], energies)]
            else:
                for i in range(len(stark_map)):
                    c = characters[i, color_mode]
                    stark_map[i][-1][1].append(energies[i])
                    if stark_map[i][-1][0] != c:
                        stark_map[i].append((c, [energies[i]]))

    GHz = 299792458 * 1e-9 * 100
    for line in stark_map:
        idx = 0
        for c, subline in line:
            plt.plot(efields[idx: idx+len(subline)], np.array(subline) * GHz, '-', color=COLOR_PALETTE[np.clip(int(c), 0, len(COLOR_PALETTE)-1), :])
            idx += len(subline) - 1

    ax.set_xlim((0, 16))
    ax.set_ylim((-4450, -4415))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $h$ (GHz)")
    ax.set_title(os.path.basename(path))
    fig.tight_layout()

    return fig, ax


def plot_starkmap_lines(path, color_mode=2):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        
        efields = np.empty(len(keys))
        stark_map = []

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efields[idx] = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            states = np.array(datagroup["States"][:,:])
            characters = np.array(datagroup["Character"])

            if idx == 0:
                stark_map = [[(c, [e])] for c, e in zip(characters[:, color_mode], energies)]
            else:
                for i in range(len(stark_map)):
                    c = characters[i, color_mode]
                    stark_map[i][-1][1].append(energies[i])
                    if stark_map[i][-1][0] != c:
                        stark_map[i].append((c, [energies[i]]))

    GHz = 299792458 * 1e-9 * 100
    for line in stark_map:
        idx = 0
        for c, subline in line:
            plt.plot(efields[idx: idx+len(subline)], np.array(subline) * GHz, '-', color=COLOR_PALETTE[np.clip(int(c), 0, len(COLOR_PALETTE)-1), :])
            idx += len(subline) - 1

    ax.set_xlim((0, 16))
    ax.set_ylim((-4450, -4415))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $h$ (GHz)")
    ax.set_title(os.path.basename(path))
    fig.tight_layout()

    return fig, ax


if __name__ == '__main__':

    dir_paths = {
        "ludwigsburg": "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/",
        "monaco": "/home/pneufeld/git/QuantumSimulation/build/apps/NO_StarkMap",
        "panama": "/mnt/Data/pneufeld/Masterarbeit/06_StarkMap/03_NO/",
    }
    dir_path = dir_paths[socket.gethostname()]
    paths = sorted(glob(f"{dir_path}/NOStarkMap*.h5"))[-1:]
    
    for path in paths:
        print(path)
        fig, ax = plot_starkmap_lines(path)

        # fig.savefig(os.path.join(dir_path, os.path.basename(path) + ".pdf"))
        # fig.savefig(os.path.join(dir_path, os.path.basename(path) + ".png"))

    plt.show()
