import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 15})

# Loading data from simulations
norm_dev_0 = pa.mat()
norm_dev_0.load("../files/deviation_0_slits.bin")

norm_dev_2 = pa.mat()
norm_dev_2.load("../files/deviation_2_slits.bin")

# Creating time array
t = np.arange(0, np.size(norm_dev_0))

# Plotting
plt.subplots_adjust(bottom = 0.138)
plt.grid()
plt.plot(t, norm_dev_0, label="No slits")
plt.plot(t, norm_dev_2, label="Double-slit")
plt.xlabel("Time steps")
plt.ylabel("$|1-|p_{ij}||$")
plt.legend()
plt.savefig("../figures/deviations.pdf")
plt.show()
