import os
import numpy as np
from statsmodels.stats.correlation_tools import corr_nearest


def get_npz_files(active_directory):
    """
    Returns npz files at path given
    :param active_directory:
    :return:
    """
    files = os.listdir(active_directory)

    npz_list = []

    for file in files:
        if ".npz" in file:
            npz_list.append(file)

    return npz_list


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

        for sub_dir in sub_dirs:

            if sub_dir in possible_sub_dirs:
                active_dir = sim_dir + temp_dir + "/" + sub_dir

                npz_files = get_npz_files(active_dir)

                for npz_file in npz_files:
                    if "dyna_storage.npz" in npz_file:
                        number_of_jumps += get_num_jumps_in_npz(active_dir + "/" + npz_file)

        print(number_of_jumps)
