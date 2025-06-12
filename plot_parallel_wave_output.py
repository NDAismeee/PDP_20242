import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Parameters (update if needed)
Ny = 100  # Number of columns in each file (assumed fixed)

# Step 1: Read all wave_rank*.dat files
file_list = sorted(glob.glob("wave_rank*.dat"))  # e.g., wave_rank0.dat, wave_rank1.dat, ...

if not file_list:
    raise FileNotFoundError("No wave_rank*.dat files found in the current directory.")

# Step 2: Load and concatenate data block-by-block
blocks = []

for filename in file_list:
    with open(filename, 'r') as f:
        lines = f.readlines()

    row_block = []
    for line in lines:
        if line.strip() == "":
            if row_block:
                blocks.append(row_block)
                row_block = []
        else:
            row_block.append(list(map(float, line.strip().split())))
    if row_block:
        blocks.append(row_block)

# Flatten block list into full matrix
flat_data = [row for block in blocks for row in block]

# Convert to structured arrays
flat_data = np.array(flat_data)
x = flat_data[:, 0]
y = flat_data[:, 1]
z = flat_data[:, 2]

# Determine grid size from x/y coordinates
unique_x = np.unique(x)
unique_y = np.unique(y)
Nx = len(unique_x)
Ny = len(unique_y)

# Reshape into 2D arrays
X = x.reshape(Nx, Ny)
Y = y.reshape(Nx, Ny)
Z = z.reshape(Nx, Ny)

# Step 3: Plot surface
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_title('2D Wave Equation Final State (Parallel Output)')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,y)')

plt.tight_layout()
plt.show()
