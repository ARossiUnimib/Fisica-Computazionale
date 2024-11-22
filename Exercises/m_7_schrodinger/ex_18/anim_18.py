import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
from matplotlib.colors import Normalize

# Load the .dat file without headers (space-delimited)
data = pd.read_csv('ex_18.dat', sep=r'\s+', header=None)

# Assign column names based on the structure: x, y, z_1, z_2, ..., z_n
n_z = data.shape[1] - 2  # Calculate number of z columns (all except x and y)
column_names = ['x', 'y'] + [f'z_{i+1}' for i in range(n_z)]  # Dynamically assign z columns
data.columns = column_names  # Assign these names to the DataFrame

# Extract x, y coordinates and z values
x = data['x'].values
y = data['y'].values
z_columns = [col for col in data.columns if col.startswith('z_')]  # Collect all z_n columns
z_values = data[z_columns].values.T  # Transpose to make frames (n) along rows

# Normalize the color scale across all frames to ensure the colors remain consistent
vmin = np.min(z_values)
vmax = np.max(z_values)
norm = Normalize(vmin=vmin, vmax=vmax)

# Create a grid for the heatmap (this assumes you want a regular grid)
x_grid, y_grid = np.meshgrid(np.unique(x), np.unique(y))

# Create a figure and axis for plotting
fig, ax = plt.subplots()

# Create an initial pcolormesh plot (this assumes the grid is structured well)
pcm = ax.pcolormesh(x_grid, y_grid, np.reshape(z_values[0], (len(np.unique(y)), len(np.unique(x)))), cmap='inferno', norm=norm)

# Add color bar
plt.colorbar(pcm, ax=ax, label='Heatmap Value')

# Set axis labels and title
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_title('Heatmap Animation')

# Update function for animation
def update(frame):
    # Update the pcolormesh with the new frame's z values
    pcm.set_array(np.reshape(z_values[frame], (len(np.unique(y)), len(np.unique(x)))).ravel())  # Remove the extra parenthesis
    ax.set_title(f"Heatmap Animation (Frame {frame + 1})")  # Update title with frame number
    return pcm,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(z_values), interval=200, blit=True)

# Save the animation as a file (optional)
ani.save('heatmap_animation.mp4', writer='ffmpeg', fps=10)

# Show the animation
plt.show()
