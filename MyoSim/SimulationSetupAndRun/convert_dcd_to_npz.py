import os
import sys

import MDAnalysis
import numpy as np


def read_pocket_list(res_list, dcd_univ):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd_univ.select_atoms(res).ix[0])
    return np.array(pocket_indices)


def get_files(datadir_, check_files=False):
    dyna_files_ = []
    psf_files_ = []
    npz_files_ = []
    files = os.listdir(datadir_)
    for file in files:
        if "dyna" in file:
            if "dcd" in file:
                dyna_files_.append(datadir_ + "/" + file)
                if debug:
                    print(os.path.isfile(datadir_ + "/" + file))
            elif ".npz" in str(file):
                npz_str = datadir_ + "/" + file
                npz_str = npz_str.replace(".npz", ".dcd")
                npz_files_.append(npz_str)
        elif ".psf" in file:
            psf_files_.append(datadir_ + "/" + file)
            if debug:
                print(os.path.isfile(datadir_ + "/" + file))
    if check_files:
        copy_dyna_files_ = dyna_files_.copy()
        for file in copy_dyna_files_:
            if file in npz_files_:
                dyna_files_.remove(file)
        if len(dyna_files_) > 0:
            dyna_files_.pop(-1)
        print("NPZ:: \n", npz_files_)
        print("DYNA:: \n", dyna_files_)

    return dyna_files_, psf_files_


def generate_npz(dyna_files_, path_psf_):

    if len(dyna_files_) < 1:
        return

    with open(out_filepath, "a") as file:
        file.write("\n" + "generating npz for " + str(dyna_files_[0]))
        file.close()

    print("generating npz for", dyna_files_)
    distfile = dyna_files_[0].replace(".dcd", "") + "all" + ".npz"
    dcd = MDAnalysis.Universe(path_psf_, dyna_files_)

    CO = dcd.select_atoms("resname CO")
    coords = []

    for ts in dcd.trajectory:
        coords.append(CO.positions)
    coords = np.array(coords)

    types = np.array([str(CO.types[0]), str(CO.types[1])])
    np.savez(distfile, coordinates=coords, types=types)


if __name__ == "__main__":
    print("Arguments passed:", sys.argv)

    if len(sys.argv) == 2:
        temp_dirs = [sys.argv[1]]
    else:
        temp_dirs = ["1K", "5K", "10K", "50K", "100K", "300K"]
        temp_dirs = ["50K", "100K", "300K"]
    debug = False
    exec_one = False

    CO_atoms = ["resname CO and type C", "resname CO and type O"]

    out_filepath = "converterOUT.txt"
    if not os.path.exists(out_filepath):
        with open(out_filepath, "w") as file:
            pass

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
        sim_dirs = [cwd + "/WITH_CO_V2/", cwd + "/NO_CO_V2/"]
        possible_sub_dirs = np.arange(1, 51).tolist()
        for i in range(len(possible_sub_dirs)):
            possible_sub_dirs[i] = str(possible_sub_dirs[i])

        for sim_dir in sim_dirs:
            for temp_dir in temp_dirs:
                sub_dirs = sorted(os.listdir(sim_dir + temp_dir))
                print(temp_dir)
                for sub_dir in sub_dirs:
                    if sub_dir in possible_sub_dirs:
                        active_dir = sim_dir + temp_dir + "/" + sub_dir
                        print("active_dir", active_dir)
                        dyna_files, psf_files = get_files(active_dir)
                        print("GET files complete")

                        psf_file = psf_files[0]

                        print("start gen of npz files")
                        generate_npz(dyna_files, psf_file)
                        print("Completed npz file gen")
                        with open(out_filepath, "a") as file:
                            file.write("\n" + "Writing to tempdir" + temp_dir + " subdir" + str(sub_dir))
                            file.close()