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
COLOR_PALETTE = np.array(
    [col.to_rgba(c) for c in ("black", "blue", "red", "green", "purple", "orange", "gray")]
)
energyInvCm = 100 * 299792458 * 6.626e-34
energyGHz = 6.626e-34 * 1e9


#
# Coloring functions
#

def color_rot(datagroup, basis):
    return COLOR_PALETTE[
        np.clip(datagroup["Character"][:, 2].astype(int), 0, len(COLOR_PALETTE) - 1)
    ]

def color_n42_R5(datagroup, basis):
    ns = datagroup["Character"][:, 0].astype(int)
    Rs = datagroup["Character"][:, 2].astype(int)
    return np.array([[1.0, 1.0, 1.0, 0.5 if n == 42 and R==5 else 0.0] for n, R in zip(ns, Rs)])

def color_fcharacter(datagroup, basis):
    mask = np.logical_or(basis[:, 1] == 3, basis[:, 1] == 1)
    states = datagroup["States"][:, :]
    strength = np.clip(np.sum((states[mask, :] ** 2), axis=0), 0, 1)
    return np.array([[0.0, 1.0, 0.0, x] for x in strength])


def color_fcharacter_N5(datagroup, basis):
    mask = np.logical_and(np.logical_or(basis[:, 1] == 3, basis[:, 1] == 1), basis[:, 2] == 5)
    states = datagroup["States"][:, :]
    strength = np.clip(np.sum((states[mask, :] ** 2), axis=0), 0, 1)
    return np.array([[0.0, 1.0, 0.0, np.clip(2.5*x, 0.0, 1.0)] for x in strength])

#
# Plot functions
#

def plot_stark_theory_lines(ax, path, color_func, xoff=0, yoff=0, unit=energyGHz, ecorr=1.0):
    with h5py.File(path) as file:
        keys = [k for k in file.keys() if k.isdecimal()]
        efields = np.empty(len(keys))
        stark_map = np.empty((len(file[keys[0]]["Energies"][:, 0]), len(keys)))
        colors = np.empty((*stark_map.shape, 4))
        basis = file["Basis"][:, :]

        for idx, datagroup in enumerate(tqdm.tqdm([file[k] for k in keys])):
            efields[idx] = float(datagroup.attrs["Electric_Field"]) / 100 * ecorr - xoff
            stark_map[:, idx] = np.array(datagroup["Energies"][:, 0]) / unit - yoff
            colors[:, idx, :] = color_func(datagroup, basis)
            
    # flip, TODO
    # flipoff = 3.25
    # stark_map = -stark_map + 2*np.ones_like(stark_map)*flipoff-0.7

    ymin, ymax = ax.get_ylim()
    for idx in range(len(stark_map)):

        if np.all(np.logical_or(stark_map[idx:] > ymax, stark_map[idx:] < ymin)):
            continue

        vertices = np.array([efields, stark_map[idx, :]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([vertices[:-1], vertices[1:]], axis=1)
        lc = LineCollection(segments, colors=colors[idx], linewidths=0.5)
        ax.add_collection(lc)

def stark_theory_plot(exdata, cmap_func, simdata_path, plot_path, yoff, cfunc, ecorr=1.0):
    fig, ax = stark_plot_helper(*exdata, cmap_func)

    plot_stark_theory_lines(ax, simdata_path, cfunc, -0.7, yoff - 2.4, ecorr=ecorr)

    fig.tight_layout()
    fig.savefig(plot_path)


def stark_theory_plot_standalone(simdata_path, plot_path, cfunc, ymin=-66.5, ymax=-61.5, xmin=0.0, xmax=25.0):
    fig = plt.figure(figsize=plotlib.set_size(aspect=0.484, fraction=1.0))
    ax = fig.add_subplot()

    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))

    plot_stark_theory_lines(ax, simdata_path, cfunc, unit=energyInvCm)

    ax.set_xlabel(r"Electric field (V/cm)")
    ax.set_ylabel(r"Energy/hc (cm${}^-1$)")

    fig.tight_layout()
    fig.savefig(plot_path)


#
# MAIN Functions
#

def main_overlap(data_folder, plot_folder):
    
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

    data = {
        "n42R5_1": ["NOStarkMap_20221111-121150_calcc.h5", -77.303],
        "n42R5_2": ["NOStarkMap_20221111-200322_calcc.h5", -77.303],
        "n42R5_3": ["NOStarkMap_20221117-143542_calcc.h5", -78.268],
        "n42R3": ["NOStarkMap_20221114-114554_calcc.h5", -1149.86],
        "n42R5_4": ["NOStarkMap_20221125-190505_calcc.h5", -78.268]
    }["n42R5_4"]
    filename, yoff = data

    data_path = os.path.join(data_folder, filename)
    plot_path = os.path.join(plot_folder, filename + ".pdf")

    # stark_theory_plot(meas_data, func, data_path, plot_path, yoff, color_rot)
    stark_theory_plot(meas_data, func, data_path, plot_path, yoff, color_fcharacter_N5, ecorr=4.5)

def main_hogan_test(data_folder, plot_folder):
    # name = "NOStarkMap_20221117-100108_calcc.h5"
    # name = "NOStarkMap_20221117-143542_calcc.h5"
    name = "NOStarkMap_20230101-235425_panama.h5"
    data_path = os.path.join(data_folder, name)
    plot_path = os.path.join(plot_folder, name + "_theory.pdf")
    stark_theory_plot_standalone(data_path, plot_path, color_rot, ymin=-66.5, ymax=-61.5, xmax=25)

if __name__ == "__main__":

    plotlib.plot_setup(no_pgf=True)

    app_folder = get_default_app_folder("NO_StarkMap")
    data_folder = os.path.join(app_folder, "Data")
    plot_folder = os.path.join(app_folder, "Plots")

    os.makedirs(data_folder, exist_ok=True)
    os.makedirs(plot_folder, exist_ok=True)

    # main_overlap(data_folder, plot_folder)
    main_hogan_test(data_folder, plot_folder)
    