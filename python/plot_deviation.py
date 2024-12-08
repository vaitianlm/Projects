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
plt.ylabel("Deviation from 1")
plt.show()

# print(pa.size(S)) #199x199x2

# S_py = np.array(S)
# S_py_magn = np.abs(S_py)

# print(np.shape(S_py))


# # Plotting
# plt.figure(figsize=(8, 6))
# im = plt.imshow(S_py_magn[0,:,:]**2, cmap='viridis', origin='lower', aspect='auto')
# plt.colorbar(im, label='Magnitude')  # Change label based on what you plot
# plt.title('First Slice of Quantum Box - Magnitude')
# plt.xlabel('X-axis (Columns)')
# plt.ylabel('Y-axis (Rows)')
# plt.tight_layout()
# plt.show()