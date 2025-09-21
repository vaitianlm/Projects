import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyarma as pa
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

z_data_list = pa.cx_cube()
z_data_list.load("../files/double_slit.bin")
z_data_list = pa.abs(z_data_list)
z_data_list = np.array(z_data_list)

# Some settings
fontsize = 15
t_min = 0
dt = 2.5e-5
x_min, x_max = 0, 1
y_min, y_max = 0, 1


# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalisation according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Plot the first frame
img = ax.imshow(z_data_list[80], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("$x$", fontsize=fontsize)
plt.ylabel("$y$", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("|u(x,y,t)|", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)


# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

writer = PillowWriter(fps=15)  # Adjust FPS as needed
# anim.save('../animation/animation.gif', writer=writer)

# Run the animation!
plt.show()

# Save the animation
# writer = FFMpegWriter(fps=15, bitrate=10000)
# anim.save('animation.mp4', writer=writer)

# anim.save('../animation/animation.mp4', writer="ffmpeg", bitrate=10000, fps=15)  # The fps (frames per second) sets the animation speed