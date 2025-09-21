import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
plt.rcParams.update({"font.size": 15})

# Loading data from simulation
u_arma = pa.cx_cube()
u_arma.load("../files/colourmap_2_slits.bin")

# Prints dimensions
print(pa.size(u_arma)) #199x199x81

# Converting from arma cube to numpy array
u = np.array(u_arma)
u_magn = np.abs(u) # Magnitude
u_Re = np.real(u)  # Real part
u_Im = np.imag(u)  # Imaginary part

# Prints dimensions (Note how it's inverted)
print(np.shape(u)) #(81, 199, 199)

dt = 2.5e-5 # Simulation time step for picking out slices

t0 = int(0)          # Slice index at t = 0.0
t001 = int(0.001/dt) # Slice index at t = 0.001
t002 = int(0.002/dt) # Slice index at t = 0.002

# Plotting

fig = plt.figure(figsize=(6, 12))
gs = gs.GridSpec(nrows=3, ncols=2, width_ratios=[1, 0.05], right=0.8, hspace=0.12)

x_min, x_max = 0, 1
y_min, y_max = 0, 1

times = [0.0, 0.001, 0.002]
indices = [t0, t001, t002]

for i, (idx, t_val) in enumerate(zip(indices, times)):
    # Computing norm
    norm = plt.cm.colors.Normalize(vmin=np.min(u_magn[idx,:,:]), vmax=np.max(u_magn[idx,:,:])) # Choose u_magn, u_Re or u_Im
    
    # Adding the main axes
    ax = fig.add_subplot(gs[i, 0])
    im = ax.imshow(u_magn[idx,:,:], extent=[x_min, x_max, y_min, y_max],   # Choose u_magn, u_Re or u_Im
                   cmap="viridis", origin="lower", aspect="auto", norm=norm)
    
    # Setting labels and time text
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.text(0.95, 0.95, f"t = {t_val:.3e}", color="white",
            ha="right", va="top", transform=ax.transAxes)

    # Adding the colorbar in the corresponding axes defined by GridSpec
    cax = fig.add_subplot(gs[i, 1])
    cbar = fig.colorbar(im, cax=cax)

    # cbar.set_label("$|u(x,y,t)|$")
    # cbar.set_label("Re$(u(x,y,t))$")
    # cbar.set_label("Im$(u(x,y,t))$")

# plt.savefig("../figures/colourmaps_magn.pdf")
# plt.savefig("../figures/colourmaps_Re.pdf")
# plt.savefig("../figures/colourmaps_Im.pdf")

plt.show()
