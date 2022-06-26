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

    # filename = "NOStarkMap_20220428-140003_calcc"
    # dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/"

    filename = "NOStarkMap_20220619-143956_monaco"
    dir_path = "/home/pneufeld/git/QuantumSimulation/build/apps/NO_StarkMap"
    
    path = os.path.join(dir_path, filename + ".h5")

    paths = sorted(glob(f"{dir_path}/NOStarkMap*.h5"))[-1:]

    for path in paths:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_title(path.split("/")[-1])

        color_palette = np.array([col.to_rgb(c) for c in ("black", "blue", "red", "green", "purple", "orange", "grey")])
              
        prev_states = None

        with h5py.File(path) as file:

            keys = [k for k in file.keys() if k.isdecimal() and k[-1] == '0']
            
            efields = np.empty(len(keys))
            stark_map = np.empty((file[keys[0]]["States"].shape[0], len(keys)))
            colors = color_palette[np.clip(file[keys[0]]["Character"][:, 2].astype(int), 0, len(color_palette)-1), :]

            for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
                # load data from file
                efields[idx] = float(datagroup.attrs["Electric_Field"])
                eField = np.ones(stark_map.shape[0])*efields[idx]

                energies = np.array(datagroup["Energies"][:, 0])
                states = np.array(datagroup["States"])
                character = np.array(datagroup["Character"])

                # reorder basis states to match order of previous configuration
                if prev_states is not None:

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

                    # order = list(range(states.shape[0]))
                    # order = np.max(np.abs(np.matmul(np.conj(np.transpose(prev_states)), states)), axis=0)
                    # order = np.empty(energies.shape, dtype=int)
                    # for i in range(order.shape[0]):
                    #     test = np.sum(np.abs(character - prev_states[i, :]), axis=1)
                    #     order[i] = np.argmin(test)

                    test = len(set(order))

                    energies = energies[order]
                    states = states[order, :][:, order]
                    character = character[order, :]
                    stark_map[:, idx] = energies
                # else:
                #     # prev_states = character
                #     prev_states = states
                
                stark_map[:, idx] = energies

                # colors = color_palette[np.mod(character[:, 2].astype(int), len(color_palette)), :]
                # colors = color_palette[np.clip(character[:, 3].astype(int), 0, len(color_palette)-1), :]
                
                # ax.scatter(eField, energies, 0.5, colors)
                # ax.scatter(eField, energies, 0.75, colors)

                prev_states = states
                prev_chars = character
 
        for i in range(stark_map.shape[0]):
            plt.plot(efields, stark_map[i, :], '-', color=colors[i, :])

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
