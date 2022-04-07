# Philipp Neufeld, 2021-2022

import matplotlib.pyplot as plt

# This is a analysis of the runtime the calcc 
# simulation computer needs to calculate some
# example of a NO stark map with respect to the number 
# of concurrent calculation threads used.

# It is apparent that the best runtime is achieved 
# at approx. half the number of calculation threads
# than the number of logical cores on the system (64).

# Theory:
# This might be due to the fact that the system actually has 
# 32 physical cores with hyperthreading enabled. The 2 threads
# that run on a single core share their memory and thus have
# to take turns accessing it. When processing large amounts 
# of data at a time the whole L3-cache is flushed and the relevant
# data is loaded from main memory. When the other thread now needs to
# access its data, the data that lies in the cache is useless to it
# since it belongs to the other thread thus the memory is discarded and
# loaded from main memory. This now means that the memory the first thread
# loaded in the beginning has to be loaded again (and again...).
# see "cache thrashing"

t_run = {
    5: 112.3, 
    6: 93.4, 
    8: 70.9, 
    10: 57.5, 
    15: 40.1, 
    20: 35.2, 
    25: 34.4, 
    30: 38.2, 
    35: 54.3,
    37: 78.13,
    40: 113.2,
    45: 119.76,
}

if __name__ == '__main__':

    t_run = {k: t_run[k] for k in sorted(t_run)}

    plt.xlabel("Cores")
    plt.ylabel("Runtime [s]")
    # plt.loglog()
    plt.plot(t_run.keys(), t_run.values(), 'o')
    plt.show()
