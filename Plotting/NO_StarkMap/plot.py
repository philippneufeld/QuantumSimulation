# Philipp Neufeld, 2021-2022

from typing import Any, Tuple

import os
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col


class NOStarkMap:

    def __init__(self, path) -> None:
        self.path = path

    def generate_R_matrix(self, basis: np.ndarray) -> np.ndarray:
        """Generate a matrix that shows which basis state belongs to which R manifold

        Args:
            basis (np.ndarray): basis states

        Returns:
            np.ndarray: Generated matrix
        """
        maxR = np.max([b[2] for b in basis])
        rot_matrix = np.zeros((int(maxR)+1, len(basis)))
        for i in range(len(basis)):
            rot_matrix[int(basis[i, 2]), i] = 1
        return rot_matrix

    def process_datagroup_R(self, datagroup: Any, rot_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Load the slice from 

        Args:
            datagroup (np.ndarray): HDF datagroup for the desired electric field
            rot_matrix (np.ndarray): Matrix that specifies which basis belongs to which R manifold

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: electric fields, energies, colors
        """
        energies = np.array(datagroup["Energies"][::])
        states = np.array(datagroup["States"][::])
        eField = np.ones_like(energies) * float(datagroup.attrs["Electric_Field"])

        # select color based on the R manifold
        colors = np.array([col.to_rgb(f"C{np.argmax(np.sum(rot_matrix * np.square(s), axis=1))}") for s in states.T])

        # match states to previous states
        # max_overlap = np.argmax(np.matmul(statesT, prev_states), axis=0)
        # eField = eField[max_overlap]
        # energies = energies[max_overlap]
        # colors = colors[max_overlap, :]
        # states = states[:, max_overlap]

        return eField, energies, colors

    def plot(self):

        # prepare data
        plot_data = []
        with h5py.File(self.path) as file:
            rot_matrix = self.generate_R_matrix(file["Basis"][::])
            for key in tqdm.tqdm([k for k in file.keys() if k.isdecimal() and k.endswith('1')]):
                eField, energies, colors = self.process_datagroup_R(file[key], rot_matrix)
                plot_data.append((eField, energies, colors))

        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for eField, energies, colors in plot_data:
            ax.scatter(eField, energies, 0.5, colors)

        ax.set_xlim((0, 25))
        ax.set_ylim((-66.5, -61.5))
        ax.set_xlabel("Electric field (V / cm)")
        ax.set_ylabel("Energy / $hc$ (cm${}^{-1}$)")

if __name__ == '__main__':

    filename = "NOStarkMap_20220407-162632_calcc"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/06_StarkMap/03_NO/"

    stark_map = NOStarkMap(os.path.join(dir_path, filename + ".h5"))
    stark_map.plot()

    plt.show()
