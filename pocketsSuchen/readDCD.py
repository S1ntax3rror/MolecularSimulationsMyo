import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io


def calc_pocket_pos(pocket_indices):
    positions_pocket = frame._pos[pocket_indices]
    mean_position_pocket = np.mean(positions_pocket, axis=0)
    dist_H2_pocket = np.linalg.norm(mean_position_H2 - mean_position_pocket)
    return dist_H2_pocket


def read_pocket_list(res_list):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd.select_atoms(res).ix[0])
    return np.array(pocket_indices)


split_dcd = ['.', 1]  # split a each '.' and take element [1] in filename_dcd

path_psf = '/home/kaeserj/PycharmProjects/CurveFitMorse/Data/step3_pbcsetup.psf'

datadir = '/Data/dcdFiles'

file_list_dcd = []

for i in range(1, 11):
    file_list_dcd.append('/home/kaeserj/PycharmProjects/CurveFitMorse/Data/dcdFiles/dyna.' + str(i) + '.dcd')

idcds = np.array([
    int(file_dcd.split('/')[-1].split(split_dcd[0])[split_dcd[1]])
    for file_dcd in file_list_dcd])

sort_dcd = np.argsort(idcds)
file_list_dcd = np.array(file_list_dcd)[sort_dcd]
idcds = idcds[sort_dcd]

num_timesteps = 10000

H2_atoms = ["resname H2 and resid 156 and name H1", "resname H2 and resid 156 and name H2"]

pocket_XE1 = ["resname LEU and resid 89 and name C", "resname HSD and resid 93 and name C",
              "resname LEU and resid 104 and name C", "resname PHE and resid 138 and name C",
              "resname ILE and resid 142 and name C", "resname TYR and resid 146 and name C"]

pocket_XE2 = ["resname LEU and resid 72 and name C", "resname ILE and resid 107 and name C",
              "resname SER and resid 108 and name C", "resname LEU and resid 135 and name C",
              "resname PHE and resid 138 and name C", "resname ARG and resid 139 and name C"]

pocket_XE3 = ["resname TRP and resid 7 and name C", "resname LEU and resid 76 and name C",
              "resname GLY and resid 80 and name C", "resname ALA and resid 134 and name C",
              "resname LEU and resid 137 and name C", "resname PHE and resid 138 and name C"]

pocket_XE4 = ["resname GLY and resid 25 and name C", "resname ILE and resid 28 and name C",
              "resname LEU and resid 29 and name C", "resname GLY and resid 65 and name C",
              "resname VAL and resid 68 and name C", "resname LEU and resid 72 and name C"]

# original
# pocket_XE5 = ["resname LEU and resid 29 and name C", "resname HEME and resid 1154 and name FE"]

# more accurate replacement for XE5
pocket_XE5 = ['resname VAL and resid 68 and name CA', 'resname LEU and resid 29 and name CA',
              'resname GLY and resid 65 and name CA', 'resname ILE and resid 107 and name CA',
              'resname HSD and resid 64 and name CA', 'resname LEU and resid 69 and name CA',
              'resname THR and resid 67 and name CA', 'resname ILE and resid 28 and name CA',
              'resname LEU and resid 104 and name CA', 'resname PHE and resid 43 and name CA',
              'resname GLY and resid 25 and name CA', 'resname LEU and resid 72 and name CA',
              'resname LEU and resid 32 and name CA', 'resname THR and resid 39 and name CA',
              'resname ALA and resid 71 and name CA', 'resname ILE and resid 99 and name CA',
              'resname HSD and resid 93 and name CA', 'resname LEU and resid 61 and name CA',
              'resname LYS and resid 42 and name CA', 'resname LEU and resid 89 and name CA',
              'resname TYR and resid 103 and name CA']

# got replaced by XE9
# pocket_XE6 = ['resname VAL and resid 10 and name CA', 'resname MET and resid 131 and name CA', 'resname ALA and resid 127 and name CA', 'resname GLN and resid 128 and name CA', 'resname PHE and resid 123 and name CA', 'resname VAL and resid 13 and name CA', 'resname LEU and resid 9 and name CA', 'resname ALA and resid 130 and name CA', 'resname TRP and resid 14 and name CA', 'resname LEU and resid 11 and name CA', 'resname HSD and resid 12 and name CA', 'resname LEU and resid 115 and name CA', 'resname GLY and resid 129 and name CA', 'resname ASP and resid 122 and name CA', 'resname GLY and resid 124 and name CA', 'resname ILE and resid 112 and name CA', 'resname ASN and resid 132 and name CA', 'resname ASP and resid 126 and name CA', 'resname TRP and resid 7 and name CA', 'resname GLN and resid 8 and name CA', 'resname GLU and resid 6 and name CA', 'resname HSD and resid 116 and name CA', 'resname ALA and resid 125 and name CA', 'resname ALA and resid 134 and name CA', 'resname ILE and resid 111 and name CA', 'resname PRO and resid 120 and name CA', 'resname LYS and resid 133 and name CA', 'resname HSD and resid 119 and name CA', 'resname LEU and resid 135 and name CA', 'resname LEU and resid 76 and name CA', 'resname VAL and resid 17 and name CA', 'resname ALA and resid 15 and name CA']
pocket_XE6 = ['resname TRP and resid 14 and name CA', 'resname VAL and resid 13 and name CA',
              'resname VAL and resid 10 and name CA', 'resname LEU and resid 11 and name CA',
              'resname VAL and resid 17 and name CA', 'resname GLY and resid 73 and name CA',
              'resname MET and resid 131 and name CA', 'resname HSD and resid 12 and name CA',
              'resname ALA and resid 15 and name CA', 'resname LEU and resid 69 and name CA',
              'resname LEU and resid 115 and name CA', 'resname LEU and resid 76 and name CA',
              'resname LYS and resid 16 and name CA', 'resname ILE and resid 111 and name CA',
              'resname LEU and resid 72 and name CA', 'resname GLU and resid 18 and name CA',
              'resname LEU and resid 9 and name CA', 'resname ILE and resid 112 and name CA',
              'resname PHE and resid 123 and name CA', 'resname TRP and resid 7 and name CA',
              'resname ALA and resid 134 and name CA', 'resname LEU and resid 135 and name CA',
              'resname LYS and resid 77 and name CA', 'resname HSD and resid 24 and name CA']

pocket_XE7 = ['resname ILE and resid 111 and name CA', 'resname VAL and resid 17 and name CA',
              'resname LEU and resid 69 and name CA', 'resname LEU and resid 115 and name CA',
              'resname TRP and resid 14 and name CA', 'resname HSD and resid 24 and name CA',
              'resname ILE and resid 112 and name CA', 'resname VAL and resid 114 and name CA',
              'resname GLY and resid 25 and name CA', 'resname LEU and resid 72 and name CA',
              'resname ILE and resid 28 and name CA', 'resname GLY and resid 73 and name CA',
              'resname VAL and resid 21 and name CA', 'resname ALA and resid 110 and name CA',
              'resname SER and resid 108 and name CA', 'resname VAL and resid 13 and name CA',
              'resname GLU and resid 18 and name CA', 'resname LYS and resid 16 and name CA',
              'resname MET and resid 131 and name CA', 'resname HSD and resid 113 and name CA',
              'resname VAL and resid 68 and name CA', 'resname ILE and resid 107 and name CA',
              'resname THR and resid 70 and name CA', 'resname ALA and resid 15 and name CA',
              'resname VAL and resid 10 and name CA', 'resname LEU and resid 135 and name CA',
              'resname ASP and resid 27 and name CA']

pocket_XE8 = ['resname GLY and resid 73 and name CA', 'resname LEU and resid 11 and name CA',
              'resname LEU and resid 76 and name CA', 'resname TRP and resid 14 and name CA',
              'resname LYS and resid 77 and name CA', 'resname VAL and resid 10 and name CA',
              'resname ALA and resid 74 and name CA', 'resname LEU and resid 72 and name CA',
              'resname VAL and resid 13 and name CA', 'resname ILE and resid 75 and name CA',
              'resname ALA and resid 15 and name CA', 'resname TRP and resid 7 and name CA',
              'resname HSD and resid 12 and name CA', 'resname VAL and resid 17 and name CA',
              'resname MET and resid 131 and name CA', 'resname ALA and resid 134 and name CA',
              'resname GLN and resid 8 and name CA', 'resname LEU and resid 69 and name CA',
              'resname LYS and resid 78 and name CA', 'resname LEU and resid 9 and name CA',
              'resname LYS and resid 79 and name CA', 'resname GLU and resid 18 and name CA',
              'resname THR and resid 70 and name CA', 'resname ALA and resid 71 and name CA']

pocket_XE9 = ['resname TRP and resid 14 and name CA', 'resname VAL and resid 13 and name CA',
              'resname VAL and resid 10 and name CA', 'resname MET and resid 131 and name CA',
              'resname LEU and resid 115 and name CA', 'resname PHE and resid 123 and name CA',
              'resname LEU and resid 11 and name CA', 'resname VAL and resid 17 and name CA',
              'resname GLN and resid 128 and name CA', 'resname ILE and resid 112 and name CA',
              'resname ILE and resid 111 and name CA', 'resname HSD and resid 12 and name CA',
              'resname ALA and resid 127 and name CA', 'resname ALA and resid 15 and name CA',
              'resname LYS and resid 16 and name CA', 'resname LEU and resid 9 and name CA',
              'resname ASN and resid 132 and name CA', 'resname GLY and resid 73 and name CA',
              'resname ASP and resid 122 and name CA', 'resname ALA and resid 130 and name CA',
              'resname ALA and resid 134 and name CA', 'resname LEU and resid 76 and name CA',
              'resname LEU and resid 135 and name CA', 'resname HSD and resid 119 and name CA',
              'resname VAL and resid 114 and name CA', 'resname HSD and resid 116 and name CA',
              'resname LEU and resid 69 and name CA', 'resname LEU and resid 72 and name CA',
              'resname GLU and resid 18 and name CA']

# initialize arrays which will hold the distance from H2 to the corresponding pocket
pocket_arr = np.zeros(num_timesteps)
pocket_XE1_dist_arr = np.zeros(num_timesteps)
pocket_XE2_dist_arr = np.zeros(num_timesteps)
pocket_XE3_dist_arr = np.zeros(num_timesteps)
pocket_XE4_dist_arr = np.zeros(num_timesteps)
pocket_XE5_dist_arr = np.zeros(num_timesteps)
pocket_XE6_dist_arr = np.zeros(num_timesteps)
pocket_XE7_dist_arr = np.zeros(num_timesteps)
pocket_XE8_dist_arr = np.zeros(num_timesteps)
pocket_XE9_dist_arr = np.zeros(num_timesteps)
pocket_XE10_dist_arr = np.zeros(num_timesteps)

count = -1

distfile = 'pocket_arrays.npz'

if not os.path.exists(distfile):

    # Loop over dcd files
    for ii, idcd in enumerate(file_list_dcd):

        # Open dcd file
        file_dcd = idcd
        dcd = MDAnalysis.Universe(path_psf, file_dcd)

        # Get atom indices for pocket
        XE1_pocket_indices = read_pocket_list(pocket_XE1)
        XE2_pocket_indices = read_pocket_list(pocket_XE2)
        XE3_pocket_indices = read_pocket_list(pocket_XE3)
        XE4_pocket_indices = read_pocket_list(pocket_XE4)
        XE5_pocket_indices = read_pocket_list(pocket_XE5)
        XE6_pocket_indices = read_pocket_list(pocket_XE6)
        XE7_pocket_indices = read_pocket_list(pocket_XE7)
        XE8_pocket_indices = read_pocket_list(pocket_XE8)
        XE9_pocket_indices = read_pocket_list(pocket_XE9)

        # Get atom indices for H2
        H2_arr = read_pocket_list(H2_atoms)

        # Get trajectory parameter
        # Nframes: Number of frames in dcd trajectory file
        # Nskip: equivalent to 'nsavc' in the CHARMM input file
        # dt: equivalent to 'timestep' in the CHARMM input file
        Nframes = len(dcd.trajectory)
        Nskip = int(dcd.trajectory.skip_timestep)
        dt = np.round(
            float(dcd.trajectory._ts_kwargs['dt']), decimals=8) / Nskip

        # Get atom types (When passing to ASE, ASE cannot always recognize
        # the atom element correctly, e.g. atom type 'CA' of alpha carbon
        # will be read as Calcium atom, for a fix see below)
        atoms = np.array([ai for ai in dcd._topology.names.values])
        Natoms = len(atoms)

        # Get atom masses
        masses = dcd._topology.masses.values

        # Get atom charges
        charges = dcd._topology.charges.values
        # Iterate over frames
        # Eventually, if you only want to add every nth step to the ASE trajectory
        # file, you can change from 'enumerate(dcd.trajectory)' to
        # 'enumerate(dcd.trajectory[::nth])'

        for jj, frame in enumerate(dcd.trajectory):
            count += 1

            # Get updated H2 position
            positions_H2 = frame._pos[H2_arr]
            mean_position_H2 = np.mean(positions_H2, axis=0)

            # Get updated position of Pocket
            pocket_XE1_dist_arr[count] = calc_pocket_pos(XE1_pocket_indices)
            pocket_XE2_dist_arr[count] = calc_pocket_pos(XE2_pocket_indices)
            pocket_XE3_dist_arr[count] = calc_pocket_pos(XE3_pocket_indices)
            pocket_XE4_dist_arr[count] = calc_pocket_pos(XE4_pocket_indices)
            pocket_XE5_dist_arr[count] = calc_pocket_pos(XE5_pocket_indices)
            pocket_XE6_dist_arr[count] = calc_pocket_pos(XE6_pocket_indices)
            pocket_XE7_dist_arr[count] = calc_pocket_pos(XE7_pocket_indices)
            pocket_XE8_dist_arr[count] = calc_pocket_pos(XE8_pocket_indices)
            pocket_XE9_dist_arr[count] = calc_pocket_pos(XE9_pocket_indices)

            # set current pocket to the smallest element = the pocket closest to the H2
            dist_list = [pocket_XE1_dist_arr[count], pocket_XE2_dist_arr[count], pocket_XE3_dist_arr[count],
                         pocket_XE4_dist_arr[count], pocket_XE5_dist_arr[count], pocket_XE6_dist_arr[count],
                         pocket_XE7_dist_arr[count], pocket_XE8_dist_arr[count], pocket_XE9_dist_arr[count]]

            dist_list = np.array(dist_list)
            pocket_arr[count] = np.argmin(dist_list) + 1

    np.savez(distfile, pocket_arr=pocket_arr, pocket_XE1=pocket_XE1_dist_arr,
             pocket_XE2=pocket_XE2_dist_arr, pocket_XE3=pocket_XE3_dist_arr, pocket_XE4=pocket_XE4_dist_arr,
             pocket_XE5=pocket_XE5_dist_arr, pocket_XE6=pocket_XE6_dist_arr, pocket_XE7=pocket_XE7_dist_arr,
             pocket_XE8=pocket_XE8_dist_arr, pocket_XE9=pocket_XE9_dist_arr)

else:
    # load from distfile
    print("reading data")
    data = np.load(distfile)
    pocket_arr = data['pocket_arr']
    pocket_XE1_dist_arr = data['pocket_XE1']
    pocket_XE2_dist_arr = data['pocket_XE2']
    pocket_XE3_dist_arr = data['pocket_XE3']
    pocket_XE4_dist_arr = data['pocket_XE4']
    pocket_XE5_dist_arr = data['pocket_XE5']
    pocket_XE6_dist_arr = data['pocket_XE6']
    pocket_XE7_dist_arr = data['pocket_XE7']
    pocket_XE8_dist_arr = data['pocket_XE8']
    pocket_XE9_dist_arr = data['pocket_XE9']

max_val = np.max([
    np.max(pocket_XE1_dist_arr), np.max(pocket_XE2_dist_arr), np.max(pocket_XE3_dist_arr), np.max(pocket_XE4_dist_arr),
    np.max(pocket_XE5_dist_arr), np.max(pocket_XE6_dist_arr), np.max(pocket_XE7_dist_arr)])

x = np.arange(num_timesteps)

fig, (ax1, ax2) = plt.subplots(2, 1)

# plot in which pocket the H2 is
ax1.plot(x, pocket_arr, linewidth=2.0, label='pocket state')

# plot all distances to H2 for debugging
ax2.plot(x, pocket_XE1_dist_arr, linewidth=1.0, label='dist to XE1')
ax2.plot(x, pocket_XE2_dist_arr, linewidth=1.0, label='dist to XE2')
ax2.plot(x, pocket_XE3_dist_arr, linewidth=1.0, label='dist to XE3')
ax2.plot(x, pocket_XE4_dist_arr, linewidth=1.0, label='dist to XE4')
ax2.plot(x, pocket_XE5_dist_arr, linewidth=1.0, label='dist to XE5')
ax2.plot(x, pocket_XE6_dist_arr, color='black', linewidth=2.0, label='dist to XE6')
ax2.plot(x, pocket_XE7_dist_arr, color='gray', linewidth=1.0, label='dist to XE7')
ax2.plot(x, pocket_XE8_dist_arr, color='gold', linewidth=1.0, label='dist to XE8')
ax2.plot(x, pocket_XE9_dist_arr, color='yellow', linewidth=1.0, label='dist to XE9')
ax2.legend()

plt.show()
