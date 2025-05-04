import os

import MDAnalysis
import numpy as np
from matplotlib import pyplot as plt
from multiprocessing import Pool


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


def plot_temp_dir(curr_dir):
    sub_dirs = sorted(os.listdir(curr_dir))
    print("processing directory " + curr_dir)

    possible_sub_dirs = np.arange(1, 51).tolist()

    for i in range(len(possible_sub_dirs)):
        possible_sub_dirs[i] = str(possible_sub_dirs[i])

    all_pockets = []

    for sub_dir in sub_dirs:
        if sub_dir in possible_sub_dirs:
            active_dir = curr_dir + "/" + sub_dir
            npz_files = get_npz_files(active_dir)

            for npz_file in npz_files:
                if "dyna_storage.npz" in npz_file:
                    fe_pos = get_iron_position(active_dir)
                    print(fe_pos)

                    data = np.load(active_dir + "/" + npz_file)
                    coordinates = data["coordinates"]
                    normalized_coords = coordinates - fe_pos

                    norm_list = np.linalg.norm(normalized_coords, axis=(1, 2))  # calc dist to FE

                    plt.clf()
                    plt.close()

                    # plt.ylabel("Distance to FE in Angström")
                    # plt.xlabel("Frame number")
                    # plt.hist(norm_list, bins=30)

                    # Plot histogram with improved appearance
                    plt.hist(norm_list, bins=30, edgecolor="black", alpha=0.75, color="royalblue", density=True)

                    # Labels and title
                    plt.xlabel("Distance to Fe (Å)", fontsize=12)
                    plt.ylabel("Probability", fontsize=12)
                    plt.title("Histogram of Distances to Fe", fontsize=14)

                    # Grid for better readability
                    plt.grid(axis="y", linestyle="--", alpha=0.7)

                    if live:
                        if "WITH_CO_V2" in curr_dir:
                            print("saving figure" + " pictures/CO/" + curr_dir.split("/")[-1] + active_dir.split("/")[-1] + ".png")
                            plt.savefig("pictures/CO/" + curr_dir.split("/")[-1] + active_dir.split("/")[-1] + ".png")
                        else:
                            plt.savefig("pictures/NO_CO/" + curr_dir.split("/")[-1] + active_dir.split("/")[-1] + ".png")
                            print("saving figure" + "pictures/NO_CO/" + curr_dir.split("/")[-1] + ".png")
                    else:
                        plt.show()




live = False

if live:
    cwd = os.getcwd()
else:
    cwd = os.getcwd() + "/../SimulationSetupAndRun"

sim_dirs = [cwd + "/WITH_CO_V2/", cwd + "/NO_CO_V2/"]
# temp_dirs = ["1K", "5K", "10K", "50K", "100K", "300K"]
temp_dirs = ["50K", "100K", "300K"]
if not live:
    temp_dirs = ["1K"]

total_dirs_w = []
total_dirs_n = []
total_dirs = []
for temp_dir in temp_dirs:
    total_dirs_w.append(sim_dirs[0] + temp_dir)
    total_dirs_n.append(sim_dirs[1] + temp_dir)
    total_dirs.append(sim_dirs[0] + temp_dir)
    total_dirs.append(sim_dirs[1] + temp_dir)

# with Pool(processes=len(total_dirs_w)) as pool:
#     pool.map(plot_temp_dir, total_dirs_w)
#
# with Pool(processes=len(total_dirs_n)) as pool:
#     pool.map(plot_temp_dir, total_dirs_n)


if not live:
    for _dir in total_dirs:
        plot_temp_dir(_dir)
else:
    with Pool(processes=len(total_dirs)) as pool:
        pool.map(plot_temp_dir, total_dirs)
