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


def plot_starkmap(path, naive=False):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    prev_states = None

    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        
        efields = np.empty(len(keys))
        stark_map = np.empty((file[keys[0]]["States"].shape[0], len(keys)))
        colors = COLOR_PALETTE[np.clip(file[keys[0]]["Character"][:, 2].astype(int), 0, len(COLOR_PALETTE)-1), :]

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            # load data from file
            efields[idx] = float(datagroup.attrs["Electric_Field"])
            energies = np.array(datagroup["Energies"][:, 0])
            states = np.array(datagroup["States"]) if not naive else None
            character = np.array(datagroup["Character"])

            # reorder basis states to match order of previous configuration
            if prev_states is not None and not naive:
                av_states = list(range(states.shape[0]))
                order = np.empty(energies.shape, dtype=int)
                for i in range(states.shape[0]):
                    s = states[av_states, i]
                    prev = prev_states[av_states, :][:, av_states]
                    chr_slice = [0, 1, 2]
                    state_overlap = np.abs(np.matmul(np.transpose(prev), s))
                    chr_oberlap = np.sum(np.abs(character[av_states, :][:, chr_slice] - prev_chars[i, :][chr_slice]), axis=1)
                    match = np.argmin(chr_oberlap - 0.25*state_overlap)
                    order[i] = av_states[match]
                    av_states.remove(av_states[match])

                # av_states = list(range(states.shape[0]))
                # order = np.empty(energies.shape, dtype=int)
                # for i in range(states.shape[0]):
                #     s = prev_states[av_states, i]
                #     prev = states[av_states, av_states]
                #     match = np.argmax(np.abs(np.matmul(np.transpose(prev), s)))
                #     order[i] = av_states[match]
                #     av_states.remove(av_states[match])

                assert(len(set(order)) == states.shape[0])

                energies = energies[order]
                states = states[order, :][:, order]
                character = character[order, :]
                stark_map[:, idx] = energies
            
            stark_map[:, idx] = energies

            prev_states = states
            prev_chars = character

    for i in range(stark_map.shape[0]):
        plt.plot(efields, stark_map[i, :], '-', color=colors[i, :])

    ax.set_xlim((0, 25))
    ax.set_ylim((-66.5, -61.5))
    ax.set_xlabel("Electric field (V / cm)")
    ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

    return fig, ax

if __name__ == '__main__':

    dir_paths = {
        "ludwigsburg": "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/",
        "monaco": "/home/pneufeld/git/QuantumSimulation/build/apps/NO_StarkMap",
        "panama": "/home/pneufeld/git/QuantumSimulation/build/apps/NO_StarkMap",
    }
    dir_path = dir_paths[socket.gethostname()]
    paths = sorted(glob(f"{dir_path}/NOStarkMap*.h5"))[-1:]

    for path in paths:
        fig, ax = plot_starkmap_scatter(path)
        # fig, ax = plot_starkmap(path, naive=True)
        
        fig.tight_layout()
        # fig.savefig(os.path.join(dir_path, os.path.basename(path) + ".pdf"))
        # fig.savefig(os.path.join(dir_path, os.path.basename(path) + ".png"))

    plt.show()
