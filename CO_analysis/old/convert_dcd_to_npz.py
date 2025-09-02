import os
import sys

import MDAnalysis
import numpy as np
import datetime
from statsmodels.imputation.ros import impute_ros


def read_pocket_list(res_list, dcd_univ):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd_univ.select_atoms(res).ix[0])
    return np.array(pocket_indices)


def get_files(datadir_, check_files=True):
    dyna_files_ = []
    psf_files_ = []
    npz_files_ = []
    files = os.listdir(datadir_)
    for file in files:
        if "dyna" in file:
            if "dcd" in file:
                dyna_files_.append(datadir_ + "/" + file)
                if debug and debug_mode == 2:
                    print("dcd file exists: ", os.path.isfile(datadir_ + "/" + file))
            elif ".npz" in str(file):
                npz_str = datadir_ + "/" + file
                npz_str = npz_str.replace(".npz", ".dcd")
                npz_files_.append(npz_str)
        elif ".psf" in file:
            psf_files_.append(datadir_ + "/" + file)
            if debug and debug_mode == 2:
                print("PSF file exists:", os.path.isfile(datadir_ + "/" + file))

    if debug and debug_mode == 2:
        print("-----------------before------------------")
        print("NPZ:: \n", npz_files_)
        print("DYNA:: \n", dyna_files_)

    if check_files:
        copy_dyna_files_ = dyna_files_.copy()
        for file in copy_dyna_files_:
            if file in npz_files_:
                dyna_files_.remove(file)
        if len(dyna_files_) > 0:
            dyna_files_ = sorted(dyna_files_)
            dyna_files_.pop(-1)

        if debug and debug_mode == 2:
            print("-----------------after------------------")
            print("NPZ:: \n", npz_files_)
            print("DYNA:: \n", dyna_files_)

    return dyna_files_, psf_files_


def generate_npz(dyna_files_, path_psf_):
    if len(dyna_files_) == 0:
        return

    indcs = []

    for ff in dyna_files_:
        idx = int(ff.split("dyna")[-1].split(".")[0])
        indcs.append(idx)
    isorted = np.argsort(indcs)
    dyna_files_ = np.array(dyna_files_)[isorted]
    dyna_files_ = dyna_files_.tolist()

    with open(out_filepath, "a") as file:
        file.write("\n" + "generating npz for " + str(dyna_files_[0]) + ":::: time is: " + str(datetime.datetime.now()))
        file.close()
    if debug:
        print("generating npz for", str(dyna_files_[0]))
    distfile = path_psf_[:-18] + "dyna_storage"
    distfile = distfile + ".npz"

    dcd = MDAnalysis.Universe(path_psf_, dyna_files_)

    coords = np.zeros((len(dcd.trajectory), 2, 3))

    CO = dcd.select_atoms("resname CO")
    CO_coords = dcd.select_atoms("resname CO").ix[0]
    for i, ts in enumerate(dcd.trajectory):
        coords[i] = ts._pos[CO_coords]
        if not i%1000:
            with open(out_filepath, "a") as file:
                file.write(f"Step {i:8d}/{len(dcd.trajectory):8d}, time: " + str(datetime.datetime.now()))
                file.close()
            print(f"Step {i:8d}/{len(dcd.trajectory):8d}, time: ", datetime.datetime.now())
    coords = np.array(coords)

    types = np.array([str(CO.types[0]), str(CO.types[1])])
    np.savez(distfile, coordinates=coords, types=types)


if __name__ == "__main__":
    print("Arguments passed:", sys.argv)

    if len(sys.argv) == 2:
        temp_dirs = [sys.argv[1]]
        out_filepath = "converterOUT" + str(sys.argv[1]) + ".txt"
    else:
        temp_dirs = ["1K", "5K", "10K", "50K", "100K", "300K"]
        #temp_dirs = ["50K", "100K", "300K"]
        out_filepath = "converterOUT.txt"

    debug = False
    debug_mode = 2
    exec_one = False
    test_1 = False
    check_for_npz_already_existing = False

    CO_atoms = ["resname CO and type C", "resname CO and type O"]

    with open(out_filepath, "w") as file:
        file.write("")

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
        if debug:
            print("cwd:", cwd)

        if test_1:
            sim_dirs = [cwd + "/NO_CO_V2/"]
            possible_sub_dirs = "1"
        else:
            sim_dirs = [cwd + "/WITH_CO_V2/", cwd + "/NO_CO_V2/"]
            possible_sub_dirs = np.arange(1, 51).tolist()

            for i in range(len(possible_sub_dirs)):
                possible_sub_dirs[i] = str(possible_sub_dirs[i])

        for sim_dir in sim_dirs:
            for temp_dir in temp_dirs:
                sub_dirs = sorted(os.listdir(sim_dir + temp_dir))

                for sub_dir in sub_dirs:
                    if sub_dir in possible_sub_dirs:
                        active_dir = sim_dir + temp_dir + "/" + sub_dir

                        if debug:
                            print("active_dir", active_dir)

                        dyna_files, psf_files = get_files(active_dir, check_for_npz_already_existing)

                        if debug:
                            print("GET filses complete")

                        psf_file = psf_files[-1]

                        if debug:
                            print("start gen of npz files")

                        generate_npz(dyna_files, psf_file)

                        if debug:
                            print("Completed npz file gen")

                        with open(out_filepath, "a") as file:
                            file.write("\n" + "Writing to tempdir" + temp_dir + " subdir" + str(sub_dir))
                            file.close()
