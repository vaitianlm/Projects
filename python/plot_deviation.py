import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

norm_dev = pa.mat()
norm_dev.load("files/deviation_0slit.bin")

t = np.arange(0, np.size(norm_dev))

plt.rcParams.update({"font.size": 16})
plt.subplots_adjust( bottom = 0.138)
plt.grid()
plt.plot(t, norm_dev)

plt.xlabel("Time steps")
plt.ylabel("Deviation from 1.0")
plt.savefig("figures/deviation_noslit.pdf")
plt.show()
