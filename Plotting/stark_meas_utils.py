
import tqdm
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from fitfun import poly1, poly2

# constants
FSR835 = 900.474e6
ELECTRODE_DISTANCE = 2.1  # cm
TIA_RESISTANCE = 470  # MOhm

def get_freq_axis(t, cavity, fsr, conv_factor, height=None, poly=None):
    if height is None:
        height = np.max(cavity)/2
    dist = int(cavity.size / 100)
    peaks, _ = find_peaks(cavity, height=height, distance=dist)

    # plt.plot(t, cavity)
    # plt.plot(t[peaks], cavity[peaks], 'x')
    # plt.show()

    if len(peaks) == 0:
        return None

    lp = peaks.size
    if poly is None:
        if lp > 2:
            poly = poly2
        elif lp == 2:
            poly = poly1
    freq = np.arange(len(peaks)) * fsr * conv_factor
    fpopt, _ = curve_fit(poly, t[peaks], freq)
    freq = poly(t, *fpopt)
    return freq

def normalize_arr(arr):
    return (arr - np.min(arr))/(np.max(arr)-np.min(arr))

def _get_h5data_generator(fname, tmin, tmax, notimecut, N, stride, nomean=False):
    with h5py.File(fname, "r") as h5data:
        frequency_trigger_factor = (1)  # fallende Flanke -> scan richtung kleinerer Frequenz
        tia_question_asked = False
        freqs = None
        keys = list(sorted(h5data.keys())[::stride])
        for key in tqdm.tqdm(keys):
            measurement = h5data[key]
            try:
                sweeps = measurement.attrs["sweeps"]
            except KeyError:
                sweeps = 1
            pdata = measurement["Pressure"]
            edata = measurement["Electronics"]
            voltage = edata.attrs["electrode_voltage"]

            tia_res = edata.attrs["tia_resistor"]
            if tia_res != TIA_RESISTANCE:
                if not tia_question_asked:
                    input("Will use 470MOhm resistance - continue with enter! ")
                    tia_question_asked = True
                tia_res = TIA_RESISTANCE

            current_factor = 1e9 / (tia_res * 1e6)
            pressure = (pdata.attrs["before"] + pdata.attrs["after"]) / 2

            sdata = measurement["ScopeData"]

            t = sdata["time"][:][::N]
            acsignal = sdata["ch1"][:][::N]

            ## MEAN of sweeps
            if sweeps > 1 and not nomean:
                temp_mean = np.empty((sweeps, acsignal.size))
                temp_mean[0] = acsignal
                for csweep in range(1, sweeps):
                    temp_mean[csweep] = measurement[f"ScopeData{csweep:09d}"]["ch1"][:][::N]
                acsignal = np.mean(temp_mean, axis=0)
            ## MEAN end

            trange = np.nonzero((t > tmin) & (t < tmax))
            if not notimecut:
                t = t[trange]
                acsignal = acsignal[trange]
            try:
                cavity = sdata["ch2"][:][::N]
                trigger = sdata["ch3"][:][::N]
                if not notimecut:
                    cavity = cavity[trange]
                    trigger = trigger[trange]
                freqs = get_freq_axis(t, cavity, FSR835, 1e-9)
                trig_min = np.argmin(trigger)
                trig_max = np.argmax(trigger)
                if trig_min > trig_max:
                    frequency_trigger_factor = -1
                elif trig_max > trig_min:
                    frequency_trigger_factor = 1
                freqs *= frequency_trigger_factor
                if freqs is None:
                    raise RuntimeError("No Freq, no fun")
            except Exception as e:
                ...

            yield {
                "t": t,
                "acsignal": acsignal,
                "cavity": cavity,
                "trigger": trigger,
                "voltage": voltage,
                "pressure": pressure,
                "current_factor": current_factor,
                "freqs": freqs,
            }


def _get_h5data_min_max(f, N, stride, nomean=False):
    gen = _get_h5data_generator(f, -np.inf, np.inf, True, N, stride, nomean=nomean)
    d = next(gen)
    t = d["t"]
    sig = d["acsignal"]
    sig /= np.max(sig) - np.min(sig)
    cav = d["cavity"]
    cav /= np.max(cav) - np.min(cav)
    cav += 1.5
    tri = d["trigger"]
    tri /= np.max(tri) - np.min(tri)
    tri += 3

    plt.plot(t, sig)
    plt.plot(t, cav)
    plt.plot(t, tri)
    plt.title("Click for tmin and tmax")
    a = plt.ginput(2)
    plt.close()
    gen.close()

    if a[0][0] > a[1][0]:
        return a[1][0], a[0][0]
    return a[0][0], a[1][0]


def get_h5data(fname, reset=False, notimecut=False, N=1, stride=1, nomean=False):
    if not reset:
        with h5py.File(fname, "r") as f:
            tmin = f.attrs.get("tmin", None)
            tmax = f.attrs.get("tmax", None)
            size = len(list(f.keys())[::stride])

    if reset or tmin is None:
        tmin, tmax = _get_h5data_min_max(fname, N, stride, nomean=nomean)

        with h5py.File(fname, "a") as f:
            f.attrs["tmin"] = tmin
            f.attrs["tmax"] = tmax
            size = len(list(f.keys())[::stride])

    return _get_h5data_generator(fname, tmin, tmax, notimecut=notimecut, N=N, stride=stride, nomean=nomean), size

