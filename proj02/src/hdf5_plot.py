import numpy as np
import matplotlib.pyplot as plt
import h5py

f = h5py.File('sol.h5', 'r')

data = np.array(f['numerical_solution'])
coords = np.array(f['coordinates'])
x_coords = []
y_coords = []
temps = []

for i in range(coords.shape[0]):
	for j in range(coords.shape[1]):
		x_coords.append(coords[i][j][0])
		y_coords.append(coords[i][j][1])
		temps.append(data[i][j])

plt.figure()
plt.title("2D Solution Plot, 100x100, Gauss-Seidel")
plt.xlabel("x-axis")
plt.ylabel("y-axis")
plt.tripcolor(x_coords, y_coords, temps)
plt.colorbar()
plt.savefig('sol.png')
