import os
import MDAnalysis
import numpy as np
from scipy.signal import correlate


def read_pocket_list(res_list, dcd_univ):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd_univ.select_atoms(res).ix[0])
    return np.array(pocket_indices)


live = False
debug = True

CO_atoms = ["resname CO and type C", "resname CO and type O"]

if live:
    datadir = os.getcwd()
    files = os.listdir(datadir)
else:
    datadir = "/home/kaeserj/PycharmProjects/CurveFitMorse/Data/FrozenMyoglobinPDBs"
    files = os.listdir(datadir)

dyna_files = []
psf_files = []
for file in files:
    if "dyna" in file:
        if "dcd" in file:
            dyna_files.append(datadir + "/" + file)
            if debug:
                print(os.path.isfile(datadir + "/" + file))
    if ".psf" in file:
        psf_files.append(datadir + "/" + file)
        if debug:
            print(os.path.isfile(datadir + "/" + file))

if debug:
    print(dyna_files)

path_psf = psf_files[0]

for file_dcd in dyna_files:
    distfile = file_dcd.replace(".dcd", ".npz")
    dcd = MDAnalysis.Universe(path_psf, file_dcd)

    CO = dcd.select_atoms("resname CO")
    coords = []
    for ts in dcd.trajectory:
        coords.append(CO.positions)
    coords = np.array(coords)

    types = np.array([str(CO.types[0]), str(CO.types[1])])
    np.savez(distfile, coordinates=coords, types=types)
