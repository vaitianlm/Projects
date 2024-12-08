import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

# FOR PROBLEM 9

# Load data from simulation
S_arma = pa.cx_cube()
S_arma.load("files/prob8_double.bin")

# Converting from arma cube to numpy array
S = np.array(S_arma)
P = np.abs(S)**2        # Wavefunction -> probability function


dt = 2.5e-5                 # Simulation time step for picking out slices
t_ind = int(0.002/dt)-1     # Slice indice at t = 0.002
dx = 0.005                  # Simulation stepsize in space for picking out x = 0.8
x_ind = int(0.8/dx)-1       # x-indice of x = 0.8

# P(y|x=0.8; t=0.002)
P_screen = P[t_ind, :, x_ind]/np.linalg.norm(P[t_ind, :, x_ind])
y = np.linspace(0, 1, np.size(P_screen))

# Figure to make sure it's done correctly
vis = P
vis[t_ind, :, x_ind] = 0.001
plt.figure(figsize=(8, 6))
im = plt.imshow(vis[t_ind,:,:], extent=[0, 1, 0, 1], cmap='viridis', origin='lower', aspect='auto')
plt.colorbar(im, label='Probability') 
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()

# Plotting
plt.plot(y, P_screen)
plt.show()