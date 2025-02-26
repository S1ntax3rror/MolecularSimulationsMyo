import os
import MDAnalysis
import numpy as np
from scipy.signal import correlate


def read_pocket_list(res_list, dcd_univ):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd_univ.select_atoms(res).ix[0])
    return np.array(pocket_indices)


def get_files(datadir_):
    dyna_files_ = []
    psf_files_ = []
    npz_files_ = []
    for file in files:
        if "dyna" in file:
            if "dcd" in file:
                dyna_files_.append(datadir_ + "/" + file)
                if debug:
                    print(os.path.isfile(datadir_ + "/" + file))
        elif ".psf" in file:
            psf_files_.append(datadir_ + "/" + file)
            if debug:
                print(os.path.isfile(datadir_ + "/" + file))
        elif ".npz" in file:
            npz_files_.append(datadir_ + "/" + file)
    return dyna_files_, psf_files_


def generate_npz(dyna_files_, path_psf_):
    for file_dcd in dyna_files_:
        distfile = file_dcd.replace(".dcd", ".npz")
        dcd = MDAnalysis.Universe(path_psf_, file_dcd)

        CO = dcd.select_atoms("resname CO")
        coords = []

        for ts in dcd.trajectory:
            coords.append(CO.positions)
        coords = np.array(coords)

        types = np.array([str(CO.types[0]), str(CO.types[1])])
        np.savez(distfile, coordinates=coords, types=types)


live = False
debug = True
exec_one = False

CO_atoms = ["resname CO and type C", "resname CO and type O"]

if live:
    if exec_one:
        datadir = os.getcwd()
        files = os.listdir(datadir)

        dyna_files, psf_files = get_files(datadir)

        if debug:
            print(dyna_files)

        path_psf = psf_files[0]
        generate_npz(dyna_files, path_psf)

    else:
        cwd = os.getcwd()
        sim_dirs = ["WITH_CO_V2", "NO_CO_V2"]
        temp_dirs = ["1K", "5K", "10K", "50K", "100K", "300K"]
        for sim_dir in sim_dirs:
            for temp_dir in temp_dirs:
                sub_dirs = sorted(os.listdir(temp_dir))
                for sub_dir in sub_dirs:


else:
    datadir = "/home/kaeserj/PycharmProjects/CurveFitMorse/Data/FrozenMyoglobinPDBs"
    files = os.listdir(datadir)



