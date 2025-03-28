import os

import MDAnalysis
import numpy as np
from matplotlib import pyplot as plt


def get_npz_files(active_directory):
    """
    Returns npz files at path given
    :param active_directory:
    :return:
    """
    files = os.listdir(active_directory)

    npz_list = []

    for file in files:
        if "storage.npz" in file:
            npz_list.append(file)

    return npz_list


def get_iron_position(active_dir_path):
    files = os.listdir(active_dir_path)
    psf_files = []
    dyna_files = []
    pdb_files = []

    for file in files:
        if ".dcd" in file:
            dyna_files.append(file)
        elif ".psf" in file:
            psf_files.append(file)
        elif ".pdb" in file:
            pdb_files.append(file)

    psf_file = psf_files[-1]
    if len(dyna_files) > 0:
        dcd = MDAnalysis.Universe(active_dir_path + "/" + psf_file, active_dir_path + "/" + dyna_files[0])
    elif len(pdb_files) > 0:
        dcd = MDAnalysis.Universe(psf_file, dyna_files[0])
    else:
        return None

    iron_pos = dcd.trajectory[0]._pos[dcd.select_atoms("resid 1154 and type FE").ix[0]]

    return iron_pos


def get_num_jumps_in_npz(file):
    """
    returns the number of pocket jumps that happen in a npz file
    :param file:
    :return:
    """

    data = np.load(file)
    coordinates = data["coordinates"]

    mean = np.mean(coordinates, axis=0)
    normalized_coords = coordinates - mean

    count = 0
    max_norm = 0

    for coord in normalized_coords:
        # print(coord[0][0], coord[0][1], coord[0][2])
        norm = np.linalg.norm(coord)

        if norm > 2:
            count += 1

        if max_norm < norm:
            max_norm = norm
    print("largest distance to mean position", max_norm)
    return count


live = False

if live:
    cwd = os.getcwd()
else:
    cwd = os.getcwd() + "/../SimulationSetupAndRun"

sim_dirs = [cwd + "/WITH_CO_V2/", cwd + "/NO_CO_V2/"]
temp_dirs = ["1K", "5K", "10K", "50K", "100K", "300K"]
possible_sub_dirs = np.arange(1, 51).tolist()

for i in range(len(possible_sub_dirs)):
    possible_sub_dirs[i] = str(possible_sub_dirs[i])

for sim_dir in sim_dirs:
    for temp_dir in temp_dirs:
        sub_dirs = sorted(os.listdir(sim_dir + temp_dir))
        print("processing directory " + sim_dir.replace(cwd, "") + "/" + temp_dir)

        number_of_jumps = 0
        all_pockets = None

        for sub_dir in sub_dirs:
            if sub_dir in possible_sub_dirs:
                active_dir = sim_dir + temp_dir + "/" + sub_dir
                npz_files = get_npz_files(active_dir)

                for npz_file in npz_files:
                    if "dyna_storage.npz" in npz_file:

                        fe_pos = get_iron_position(active_dir)
                        print(fe_pos)
                        number_of_jumps += get_num_jumps_in_npz(active_dir + "/" + npz_file)

                        data = np.load(active_dir + "/" + npz_file)
                        coordinates = data["coordinates"]
                        normalized_coords = coordinates - fe_pos

                        norm_list = np.linalg.norm(normalized_coords, axis=(1, 2))  # calc dist to FE

                        if all_pockets is None:
                            all_pockets = norm_list
                        else:
                            print("all pockets", all_pockets)
                            print("norm list", norm_list)
                            print(len(all_pockets))
                            all_pockets = np.concatenate((all_pockets, norm_list))
                            print(len(all_pockets))

        if not (all_pockets is None):
            plt.clf()
            plt.close()

            plt.ylabel("Distance to FE in Angström")
            plt.xlabel("Frame number")
            plt.plot(np.arange(len(all_pockets)), all_pockets, linewidth=0.05)

            if "WITH_CO_V2" in sim_dir:
                plt.savefig(cwd + "/W_CO_" + temp_dir + ".png")
            else:
                plt.savefig(cwd + "/NO_CO_" + temp_dir + ".png")
            plt.show()

        #print(number_of_jumps)
