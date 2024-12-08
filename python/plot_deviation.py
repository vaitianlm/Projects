import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

# S = pa.cx_cube()
# S.load("files/double_slit_test.bin")

norm_dev = pa.mat()
norm_dev.load("files/deviation.bin")

t = np.arange(0, np.size(norm_dev))

plt.plot(t, norm_dev)
plt.show()

print(pa.size(S)) #199x199x2

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