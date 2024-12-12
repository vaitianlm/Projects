import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 15})

# List of files and labels for each dataset
files = ["../files/detector_1_slits.bin",
         "../files/detector_2_slits.bin",
         "../files/detector_3_slits.bin"]

labels = ["Single-slit", "Double-slit", "Triple-slit"]

# Parameters
dt = 2.5e-5                 # Simulation time step for picking out slices
t_ind = int(0.002/dt)-1     # Slice indice at t = 0.002
dx = 0.005                  # Simulation stepsize in space for picking out x = 0.8
x_ind = int(0.8/dx)-1       # x-indice of x = 0.8


for file, label in zip(files, labels):
    
    # Loading data from simulation
    S_arma = pa.cx_cube()
    S_arma.load(file)

    # Converting from arma cube to numpy arrays
    S = np.array(S_arma)

    # Converting wavefunctions to probability distributions
    P = np.abs(S)**2

    # P(y|x=0.8; t=0.002)
    P_screen = P[t_ind, :, x_ind]/np.sum(P[t_ind, :, x_ind])

    # Creating y array
    y = np.linspace(0, 1, np.size(P_screen))

    # Plotting each distribution
    plt.plot(y, P_screen, label=label)

plt.subplots_adjust(left = 0.145, bottom = 0.117)
plt.grid()
plt.xlabel("$y$")
plt.ylabel("$p(y|x=0.8; t=0.002)$")
plt.legend()
# plt.savefig("../figures/detector.pdf")
plt.show()