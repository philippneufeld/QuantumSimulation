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


def plot_starkmap_lines(path):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    prev_states = None

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal() and 0 < int(k) <= 100]
        
        efields = np.empty(len(keys))
        stark_map = np.empty((file[keys[0]]["States"].shape[0], len(keys)))
        stark_map[:, :] = np.NaN
        colors = COLOR_PALETTE[np.clip(file[keys[0]]["Character"][:, 2].astype(int), 0, len(COLOR_PALETTE)-1), :]

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efields[idx] = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            states = np.array(datagroup["States"])

            # reorder basis states to match order of previous configuration
            if idx == 0:
                prev_states = states
                stark_map[:, idx] = energies
            else:
                
                overlap = np.abs(np.matmul(np.transpose(states), prev_states))
                match = np.max(overlap, axis=0)
                match_idx = np.argmax(overlap, axis=0)

                accept = match > np.sqrt(0)
                stark_map[accept, idx] = energies[match_idx[accept]]
                prev_states[:, accept] = states[:, match_idx[accept]]

                reject = np.logical_not(accept)
                stark_map[reject, idx] = np.NaN
                prev_states[:, reject] = 0

    for i in range(stark_map.shape[0]):
        plt.plot(efields, stark_map[i, :], '-', color=colors[i, :])

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

    return fig, ax

def plot_starkmap_lines2(path):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    prev_states = None

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal() and 0 < int(k) <= np.inf]
        
        efields = np.empty(len(keys))
        stark_map = []

        min_matches = np.inf

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efields[idx] = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            states = np.array(datagroup["States"])
            characters = np.array(datagroup["Character"])

            # reorder basis states to match order of previous configuration
            if idx == 0:
                prev_states = states
                stark_map = [[(c, [e])] for c, e in zip(characters[:, 2], energies)]
            else:
                
                overlap = np.abs(np.matmul(np.transpose(states), prev_states))
                match_idx = np.argmax(overlap, axis=0)

                min_matches = min(min_matches, len(set(match_idx)))

                for i, midx in enumerate(match_idx):
                    c = characters[midx, 2]
                    stark_map[i][-1][1].append(energies[midx])
                    if stark_map[i][-1][0] != c:
                        stark_map[i].append((c, [energies[midx]]))

                prev_states = states[:, match_idx]

    for line in stark_map:
        idx = 0
        for c, subline in line:
            plt.plot(efields[idx: idx+len(subline)], subline, color=COLOR_PALETTE[int(c), :])
            idx += len(subline) - 1

    print(min_matches)

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

    return fig, ax

def plot_starkmap_lines3(path):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal() and 0 <= int(k) <= 800]
        
        efields = np.empty(len(keys))
        stark_map = []

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efields[idx] = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            characters = np.array(datagroup["Character"])

            # reorder basis states to match order of previous configuration
            if idx == 0:
                stark_map = [[(c, [e])] for c, e in zip(characters[:, 2], energies)]
            else:
                for i in range(len(stark_map)):
                    c = characters[i, 2]
                    stark_map[i][-1][1].append(energies[i])
                    if stark_map[i][-1][0] != c:
                        stark_map[i].append((c, [energies[i]]))

    for line in tqdm.tqdm(stark_map):
        idx = 0
        for c, subline in line:
            plt.plot(efields[idx: idx+len(subline)], subline, '-', color=COLOR_PALETTE[int(c), :])
            idx += len(subline) - 1

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

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
        # fig, ax = plot_starkmap_scatter(path)
        # fig, ax = plot_starkmap_lines3(path)
        fig, ax = plot_starkmap_lines2(path)

        fig.tight_layout()
        # fig.savefig(os.path.join(dir_path, os.path.basename(path) + ".pdf"))
        # fig.savefig(os.path.join(dir_path, os.path.basename(path) + ".png"))

    plt.show()
