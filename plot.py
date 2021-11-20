
import numpy as np
import matplotlib.pyplot as plt

det, coeff = np.genfromtxt("build/apps/Rb87_EIT/data.txt", unpack=True)
plt.plot(det, coeff)
plt.show()
