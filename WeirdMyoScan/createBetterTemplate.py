import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io
from ase.io import read


def build_H_string(coords_to_build):
    returnval = "H  "
    #for x in range(3):
    #if coords_to_build[x] >= 0:
    returnval += "{0: 18.16}     {1: 18.16}     {2: 18.16}".format(
        *coords_to_build
        )  #str(coords_to_build[x]) + "     "
    #else:
    #    returnval += str(coords_to_build[x]) + "    "
    return returnval


# Open xyz file
path = "/WeirdMyoScan/MyoTemplate.xyz"

atoms = read(path)
nitrogen_atoms = [atom for atom in atoms if atom.symbol == 'N']
nitrogen_atoms_coords = []
for coords in nitrogen_atoms:
    nitrogen_atoms_coords.append(np.array(coords.position))
nitrogen_atoms_coords = np.array(nitrogen_atoms_coords)

print("Nitrogen atoms found:", nitrogen_atoms)
print("At Coords:", nitrogen_atoms_coords)

mean_nitrogen_atoms = np.mean(nitrogen_atoms_coords, axis=0)

print("mean:", mean_nitrogen_atoms)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

H_atoms = [atom for atom in atoms if atom.symbol == 'H']
n = len(H_atoms)
H_atoms = H_atoms[n - 2:]
H_atoms_coords = []
for coords in H_atoms:
    H_atoms_coords.append(np.array(coords.position))
H_atoms_coords = np.array(H_atoms_coords)

print("H atoms found:", H_atoms)
print("At Coords:", H_atoms_coords)

mean_H_atoms = np.mean(H_atoms_coords, axis=0)

print("mean:", mean_H_atoms)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

iron_atom = [atom for atom in atoms if atom.symbol == 'Fe']
iron_atom_coords = iron_atom[0].position

print("Iron atoms found:", iron_atom)
print("At Coords:", iron_atom_coords)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

z_vect_FE = iron_atom_coords - mean_nitrogen_atoms

print("z vector with Nitrogens", z_vect_FE)

norm = np.sqrt(np.sum(z_vect_FE ** 2))
z_vect_FE_norm = z_vect_FE / norm

print("normalized z vector", z_vect_FE_norm)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

z_vect_H = mean_H_atoms - iron_atom_coords

print("z vector with H atoms", z_vect_H)

norm = np.sqrt(np.sum(z_vect_H ** 2))
z_vect_H_norm = z_vect_H / norm

print("normalized z vector", z_vect_H_norm)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

print("correlation of z vect dependent on H vs z vect dependent on nitrogen:", np.dot(z_vect_H_norm, z_vect_FE_norm))
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

lin_independent_vector = np.array([-z_vect_H_norm[1], z_vect_H_norm[0], 0])
orthogonal_vector = np.cross(z_vect_H_norm, lin_independent_vector)
norm = np.sqrt(np.sum(orthogonal_vector ** 2))
orthogonal_vector = orthogonal_vector / norm

print("Vector orthogonal to the given vector:", orthogonal_vector)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

scan_steps_in_z_dir = 30
scan_dist_in_z_dir = 3
scan_step_size_in_z_dir = scan_dist_in_z_dir / scan_steps_in_z_dir

scan_steps_in_d_dir = 20
scan_dist_in_d_dir = 2
scan_step_size_in_d_dir = scan_dist_in_d_dir / scan_steps_in_d_dir
h2_pos_list = []

dist_h2_FE = []
dist_h2 = []

original_dist_h2 = np.linalg.norm(H_atoms_coords[1] - H_atoms_coords[0])
original_dist_h2_FE = np.linalg.norm(iron_atom_coords - (H_atoms_coords[0] + H_atoms_coords[1]) / 2)


path = "/WeirdMyoScan/distancesH2AndFeH2.txt"

with open(path, "w") as file:

    for z_dir in range(-int(scan_steps_in_z_dir//3), int((scan_steps_in_z_dir*2)//3)):
        for d_dir in range(-int(scan_steps_in_d_dir//8), int((7*scan_steps_in_d_dir)//8)):
            print(original_dist_h2_FE + z_dir * scan_step_size_in_d_dir, original_dist_h2 + d_dir * scan_step_size_in_z_dir)
            file.write(str(original_dist_h2_FE + z_dir * scan_step_size_in_d_dir) + "  " + str(original_dist_h2 + d_dir * scan_step_size_in_z_dir) + "\n")
            dist_h2_FE.append(original_dist_h2_FE + z_dir * scan_step_size_in_d_dir)
            dist_h2.append(original_dist_h2 + d_dir * scan_step_size_in_z_dir)
            h2_pos_list.append([[H_atoms_coords[
                                     0] + z_vect_H_norm * z_dir * scan_step_size_in_z_dir + orthogonal_vector * d_dir * scan_step_size_in_d_dir],
                                [H_atoms_coords[
                                     1] + z_vect_H_norm * z_dir * scan_step_size_in_z_dir - orthogonal_vector * d_dir * scan_step_size_in_d_dir]])

            # print(H_atoms_coords[0] + z_vect_H_norm * z_dir + orthogonal_vector * d_dir)
            # print(H_atoms_coords[1] + z_vect_H_norm * z_dir - orthogonal_vector * d_dir)

    h2_pos_list = np.array(h2_pos_list).reshape(-1, 2,3)
print(
    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")


template_file_path = "/WeirdMyoScan/MyoQMMMTemplate.xyz"
input_file_path = "/WeirdMyoScan/input/input_orca_scan_0000"

input_file_path_prefix = "/WeirdMyoScan/input/input_orca_scan_"

with open(template_file_path, 'r') as file:
    temp_lines = file.read()

tag1 = "replace with H1"
tag2 = "replace with H2"

for i in range(len(h2_pos_list)):
    input_file_path = input_file_path_prefix + str(i)

    lines = temp_lines[:]

    coords_H1 = h2_pos_list[i][0]
    coords_H2 = h2_pos_list[i][1]

    replacement = build_H_string(coords_H1)
    lines = lines.replace(tag1, replacement)
    replacement = build_H_string(coords_H2)
    lines = lines.replace(tag2, replacement)

    with open(input_file_path, 'w') as file:
        file.write(lines)



# import subprocess
# subprocess.run(['orca', orca_file_path])