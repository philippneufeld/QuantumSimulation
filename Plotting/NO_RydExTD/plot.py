# Philipp Neufeld, 2021-2022

from typing import Any, Tuple

import os
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt


class GasSensorTimeDependence:

    def __init__(self, path) -> None:
        self.path = path


if __name__ == '__main__':

    filename = "NORydExTD_20220427-170218_ludwigsburg"
    dir_path = "/home/PI5/pneufeld/remote_home/Masterarbeit/07_TimeDependence/01_Ion_modG"
    path = os.path.join(dir_path, filename + ".h5")

    with h5py.File(path) as file:
        for datagroup in tqdm.tqdm([file[k] for k in file.keys() if k.isdecimal()]):
            pass
