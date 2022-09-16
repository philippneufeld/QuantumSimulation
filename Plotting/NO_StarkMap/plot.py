# Philipp Neufeld, 2021-2022

import os
import socket
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.collections import LineCollection
from glob import glob

energyInvCm = 100 * 299792458 * 6.626e-34
energyGHz = 6.626e-34 * 1e9

COLOR_PALETTE = np.array([col.to_rgba(c) for c in ("black", "blue", "red", "green", "purple", "orange", "grey")])

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

def color_rot(datagroup, basis):
    return COLOR_PALETTE[np.clip(datagroup["Character"][:, 2].astype(int), 0, len(COLOR_PALETTE)-1)]

def color_fcharacter(datagroup, basis):
    mask = np.logical_or(basis[:, 1] == 3, basis[:, 1] == 1)
    states = datagroup["States"][:,:]
    strength = np.clip(np.sum((states[mask,:]**2), axis=0), 0, 1)
    return np.array([[1.0, 0.0, 0.0, x] for x in strength])

def plot_starkmap_lines(path, color_func, unit=energyInvCm, ymin=None, ymax=None):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    name = os.path.basename(path)
    date, time = [int(x) for x in name.split('_')[1].split('-')]
    
    yyyy, mm, dd = date // 10000, (date // 100) % 100, date % 100
    h, m, s = time // 10000, (time // 100) % 100, time % 100

    ax.set_title(f"Stark map {dd:02d}.{mm:02d}.{yyyy:04d} {h:02d}:{m:02d}:{s:02d}")

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        efields = np.empty(len(keys))
        stark_map = np.empty((len(file[keys[0]]["Energies"][:, 0]), len(keys)))
        colors = np.empty((*stark_map.shape, 4))
        basis = file["Basis"][:,:]

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            efields[idx] = float(datagroup.attrs["Electric_Field"]) / 100
            stark_map[:, idx] = np.array(datagroup["Energies"][:, 0]) / unit
            colors[:, idx, :] = color_func(datagroup, basis)
            
    for idx in range(len(stark_map)):
        vertices = np.array([efields, stark_map[idx,:]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([vertices[:-1], vertices[1:]], axis=1)
        lc = LineCollection(segments, colors=colors[idx], linewidths=0.5)
        ax.add_collection(lc)

    if ymin is None:
        ymin = np.min(stark_map)
    if ymax is None:
        ymax = np.max(stark_map)

    ax.set_xlabel("Electric field (Vcm^-1)")
    ax.set_ylabel("Energy")

    ax.set_xlim((np.min(efields), np.max([*efields])))
    ax.set_ylim((ymin, ymax))

    return fig, ax


if __name__ == '__main__':

    dir_paths = {
        "ludwigsburg": "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/",
        "monaco": "/home/pneufeld/git/QuantumSimulation/build/apps/NO_StarkMap",
        "panama": "/mnt/Data/pneufeld/Masterarbeit/06_StarkMap/03_NO/",
    }
    dir_path = dir_paths[socket.gethostname()]
    paths = [
        # "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/NOStarkMap_20220719-110258_ludwigsburg.h5", # no additional terms
        # "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/NOStarkMap_20220719-110542_ludwigsburg.h5", # dipole term
        *sorted(glob(f"{dir_path}/NOStarkMap*.h5"))[-1:]
    ]
    
    for path in paths:
        print(path)
        # fig, ax = plot_starkmap_lines(path, color_fcharacter, energyInvCm)
        fig, ax = plot_starkmap_lines(path, color_rot, energyGHz) # , -66.5, -61.5)

    plt.show()
