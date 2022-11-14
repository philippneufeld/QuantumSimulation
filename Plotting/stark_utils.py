import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import plotlib
from stark_meas_utils import get_h5data, ELECTRODE_DISTANCE

def sigmoid(x, a, x0):
    return 1 / (1 + np.exp(-a * (x - x0)))

def combined_sigmoids(x, a1, x01, a2, x02, w):
    return w * sigmoid(x, a1, x01) + (1 - w) * sigmoid(x, a2, x02)

def create_cmap_funcs(func, min_sig, max_sig):
    o, a = func(0), 1/(func(1)-func(0))
    stark_cmap_forward = lambda x: a*(func((x-min_sig)/(max_sig-min_sig))-o)

    x = np.linspace(min_sig, max_sig, 10000)
    fx = stark_cmap_forward(x)
    stark_cmap_backward = lambda y, x=x, fx=fx: np.interp(y, fx, x)

    return stark_cmap_forward, stark_cmap_backward

def reduce_data(x, y, N, kfactor=3):
    ks = int(kfactor*np.max([1, len(y) / N]))
    yav = np.convolve(y, np.ones(ks)/ks, mode='same')
    x0 = np.linspace(np.min(x), np.max(x), N)
    return x0, np.interp(x0, x, yav)

def load_stark_data_h5(file, N=50, start_at = None, reset=False, stride=1, nomean=False):
    data, data_size = get_h5data(file, reset=reset, N=N, stride=stride, nomean=nomean)

    # accumulate data
    voltages = []
    acsignals = []
    freqs = []
    for cnt, subdata in enumerate(data):
        freqs = subdata["freqs"]
        voltage = subdata["voltage"]
        if start_at is not None and voltage/ELECTRODE_DISTANCE < start_at:
            continue
        voltages.append(voltage)
        signal = np.array(subdata["acsignal"] * subdata["current_factor"])
        freqs, signal = reduce_data(freqs, signal, 1000)
        acsignals.append(signal)
    acsignals = np.array(acsignals)
    voltages = np.array(voltages)

    return voltages, freqs, np.transpose(acsignals)
    
def stark_plot_helper(voltages, freqs, acsignals, func=lambda x:x, aspect=plotlib.GOLDEN_RATIO, fraction=1):

    fig = plt.figure(figsize=plotlib.set_size(aspect=aspect, fraction=fraction))

    vmin, vmax = np.min(acsignals), np.max(acsignals)

    fields = np.array(voltages) / ELECTRODE_DISTANCE

    ax = fig.add_subplot(111)
    cmap = "inferno"
    forward, backward = create_cmap_funcs(func, vmin, vmax)
    cnorm = mpl.colors.FuncNorm((forward, backward), vmin=vmin, vmax=vmax)
    pcm = ax.pcolormesh(fields, freqs, acsignals, shading="auto", 
            cmap=cmap, norm=cnorm, rasterized=True)

    ax.set_xlim((np.min(fields), np.max(fields)))
    ax.set_ylim((np.min(freqs), np.max(freqs)))

    ax.set_xlabel(r"Field (V/cm)")
    ax.set_ylabel("Frequency (GHz)")

    xticks = np.floor(100*np.clip(backward(np.linspace(0,1, 6)), 0, vmax))/100
    yticks = forward(xticks)

    bar = fig.colorbar(pcm, ax=ax, ticks=xticks)
    bar.set_label("Current (nA)")

    return fig, ax

