
import numpy as np
import matplotlib.pyplot as plt

det, coeff = np.genfromtxt("build/apps/Test/data.txt", unpack=True)
plt.plot(det, coeff)
plt.show()
