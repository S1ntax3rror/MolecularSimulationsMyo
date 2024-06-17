import cubetools

import numpy as np
from matplotlib import pyplot as plt


def read_cube(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Extract number of atoms and grid dimensions
    num_atoms = int(lines[2].split()[0])
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    x_voxels = int(lines[3].split()[0])
    y_voxels = int(lines[4].split()[0])
    z_voxels = int(lines[5].split()[0])

    # Extract voxel vectors
    x_vector = np.array([float(x) for x in lines[3].split()[1:4]])
    y_vector = np.array([float(x) for x in lines[4].split()[1:4]])
    z_vector = np.array([float(x) for x in lines[5].split()[1:4]])

    # Extract atomic positions and charges
    atoms = []
    for i in range(6, 6 + num_atoms):
        atom_data = lines[i].split()
        atoms.append([int(atom_data[0]), float(atom_data[2]), float(atom_data[3]), float(atom_data[4])])

    # Extract the density data
    density_data = []
    for i in range(6 + num_atoms, len(lines)):
        density_data.extend([float(val) for val in lines[i].split()])

    # Reshape density data into the correct 3D grid
    density_data = np.array(density_data).reshape((x_voxels, y_voxels, z_voxels))

    return origin, x_vector, y_vector, z_vector, atoms, density_data


import numpy as np
from scipy.interpolate import RegularGridInterpolator


def get_plane_slice(density_data, origin, x_vector, y_vector, z_vector, point1, point2, num_points=100):
    # Calculate the normal vector of the plane
    normal = np.cross(point2 - point1, np.array([1.5, 0, 0]))  # Using an arbitrary vector for cross product

    # Calculate the grid points
    x = np.linspace(point1[0], point1[0] + (density_data.shape[0] - 1) * x_vector[0], density_data.shape[0])
    y = np.linspace(point1[1], point1[1] + (density_data.shape[1] - 1) * y_vector[1], density_data.shape[1])
    z = np.linspace(point1[2], point1[2] + (density_data.shape[2] - 1) * z_vector[2], density_data.shape[2])

    interpolator = RegularGridInterpolator((x, y, z), density_data, bounds_error=False, fill_value=None)

    # Define the plane by creating a grid of points along the plane
    plane_points = []
    for i in range(num_points):
        t = i / (num_points - 1)
        point = point1 + t * (point2 - point1)
        plane_points.append(point)

    plane_points = np.array(plane_points)

    # Interpolate the values at these points
    slice_values = interpolator(plane_points)

    # Reshape the interpolated values to form a 2D plane
    slice_2d = slice_values.reshape((num_points, num_points))

    return slice_2d


file_path = 'h2_esp.cube'
origin, x_vector, y_vector, z_vector, atoms, density_data = read_cube(file_path)

step_size = x_vector[0]
origin_x = origin[0]
n = int(abs(origin_x//step_size))

plane = density_data[n+1,:,:]

x_arr = np.arange(len(plane))*step_size*0.529
y_arr = np.arange(len(plane[0]))*y_vector[1]*0.529

x_grid, y_grid = np.meshgrid(x_arr, y_arr)


cut = 0.5
plane[plane > cut] = cut
plane[plane < 0] = 0.000001
plane * 627.5
plane = np.log10(plane)

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 28

plt.figure(figsize=(10, 8))

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

# plt.contourf(x_grid, y_grid, plane.T, 10, cmap="inferno_r")
plt.contourf(x_grid, y_grid, plane.T, 10, cmap="inferno")
plt.colorbar(label='\nlog$_1$$_0$ ESP + 1 (kcal/mol)')

plt.xlabel("x-coordinate in $\mathrm{\AA}$")
plt.ylabel("y-coordinate in $\mathrm{\AA}$")


plt.savefig(dpi=200, fname="contourH2QuadPolLog10", )

plt.show()
