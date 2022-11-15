import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.collections import LineCollection
import h5py
import tqdm

import plotlib
from directories import get_miccell_folder, get_default_app_folder
from stark_utils import load_stark_data_h5, combined_sigmoids, stark_plot_helper

# Constants
PLOT_DIR = "Plots"
COLOR_PALETTE = np.array([col.to_rgba(c) for c in ("red",)])
COLOR_PALETTE = np.array(
    [col.to_rgba(c) for c in ("blue", "green", "red", "pink", "brown")]
)
energyInvCm = 100 * 299792458 * 6.626e-34
energyGHz = 6.626e-34 * 1e9


def color_rot(datagroup, basis):
    return COLOR_PALETTE[
        np.clip(datagroup["Character"][:, 2].astype(int), 0, len(COLOR_PALETTE) - 1)
    ]


def color_n42(datagroup, basis):
    ns = datagroup["Character"][:, 0].astype(int)
    return np.array([[1.0, 0.0, 0.0, 0.5 if n == 42 else 0.0] for n in ns])


def color_fcharacter(datagroup, basis):
    mask = np.logical_or(basis[:, 1] == 3, basis[:, 1] == 1)
    states = datagroup["States"][:, :]
    strength = np.clip(np.sum((states[mask, :] ** 2), axis=0), 0, 1)
    return np.array([[1.0, 1.0, 1.0, x] for x in strength])


def plot_stark_theory_lines(ax, path, color_func, xoff=0, yoff=0, unit=energyGHz):
    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        efields = np.empty(len(keys))
        stark_map = np.empty((len(file[keys[0]]["Energies"][:, 0]), len(keys)))
        colors = np.empty((*stark_map.shape, 4))
        basis = file["Basis"][:, :]

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            efields[idx] = float(datagroup.attrs["Electric_Field"]) / 100 * 4.5 - xoff
            stark_map[:, idx] = np.array(datagroup["Energies"][:, 0]) / unit - yoff
            colors[:, idx, :] = color_func(datagroup, basis)

    ymin, ymax = ax.get_ylim()
    for idx in range(len(stark_map)):

        if np.all(np.logical_or(stark_map[idx:] > ymax, stark_map[idx:] < ymin)):
            continue

        vertices = np.array([efields, stark_map[idx, :]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([vertices[:-1], vertices[1:]], axis=1)
        lc = LineCollection(segments, colors=colors[idx], linewidths=0.5)
        ax.add_collection(lc)


def stark_theory_plot(exdata, cmap_func, simdata_path, plot_path, yoff, cfunc):
    fig, ax = stark_plot_helper(*exdata, cmap_func)

    plot_stark_theory_lines(ax, simdata_path, cfunc, -0.7, yoff - 3.4)

    fig.tight_layout()
    fig.savefig(plot_path)


if __name__ == "__main__":

    plotlib.plot_setup(no_pgf=True)

    # load measurement data
    func = lambda x: combined_sigmoids(x, 35, 0.14, 17.5, 0.55, 0.7)
    measurement_folder = os.path.join(
        get_miccell_folder(),
        "TraceGasSensing/Measurements/2022/2022-Stark-Map/Data/Voltage",
    )
    meas_path = os.path.join(
        measurement_folder, "2022-07-13_15-54-46_834.91_1.30_32.00_0.10.h5"
    )
    meas_data = load_stark_data_h5(meas_path, N=40)

    app_folder = get_default_app_folder("NO_StarkMap")
    data_folder = os.path.join(app_folder, "Data")
    plot_folder = os.path.join(app_folder, "Plots")

    os.makedirs(data_folder, exist_ok=True)
    os.makedirs(plot_folder, exist_ok=True)

    data = {
        "n42R5_1": ["NOStarkMap_20221111-121150_calcc.h5", -77.303],
        "n42R5_2": ["NOStarkMap_20221111-200322_calcc.h5", -77.303],
        "n42R3": ["NOStarkMap_20221114-114554_calcc.h5", -1149.86],
    }["n42R5_2"]
    filename, yoff = data

    data_path = os.path.join(data_folder, filename)
    plot_path = os.path.join(plot_folder, filename + ".pdf")

    stark_theory_plot(meas_data, func, data_path, plot_path, yoff, color_fcharacter)
