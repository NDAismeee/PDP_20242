import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from file
filename = "wave_output.dat"
data = np.loadtxt(filename)

# Infer grid size
# wave_output.dat is expected to be formatted in blocks with empty lines between rows
with open(filename, 'r') as f:
    lines = f.readlines()

# Count number of blocks (i.e., lines between empty lines)
rows = []
row = []
for line in lines:
    if line.strip() == '':
        if row:
            rows.append(row)
            row = []
    else:
        row.append(list(map(float, line.strip().split())))
if row:
    rows.append(row)  # last block

# Convert to numpy arrays
X = np.array([[p[0] for p in r] for r in rows])
Y = np.array([[p[1] for p in r] for r in rows])
Z = np.array([[p[2] for p in r] for r in rows])

# Plot surface
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_title('2D Wave Equation Final State')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,y)')

plt.tight_layout()
plt.show()
