import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

# FOR PROBLEM 8 (REMEMBER TO NOT WHAT WE ARE PLOTTING i.e. sqrt(p) or p)

# Load data from simulation
S_arma = pa.cx_cube()
S_arma.load("files/no_slits.bin")

# Prints dimentions
print(pa.size(S_arma)) #199x199x2

# Converting from arma cube to numpy array
S = np.array(S_arma)
S_magn = np.abs(S) # Magnitude
S_Re = np.real(S)  # Real part
S_Im = np.imag(S)  # Imaginary part


# Prints dimentions (Note how it's inverted)
print(np.shape(S))

dt = 2.5e-5 # Simulation time step for picking out slices

t001 = int(0.001/dt) # Slice indice at t = 0.001
t002 = int(0.002/dt) # Slice indice at t = 0.002


# Plotting
plt.figure(figsize=(8, 6))
im = plt.imshow(S_magn[t002,:,:]**2, cmap='viridis', origin='lower', aspect='auto')
plt.colorbar(im, label='Magnitude')  # Change label based on what you plot
plt.title('First Slice of Quantum Box - Magnitude')
plt.xlabel('X-axis (Columns)')
plt.ylabel('Y-axis (Rows)')
plt.tight_layout()
plt.show()