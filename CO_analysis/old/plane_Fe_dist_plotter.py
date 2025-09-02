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
    psf_file = None
    dyna_file = None

    for file in files:  # get the first dyna file
        if ".dcd" in file:
            dyna_file = os.path.join(active_dir_path, file)
            break
    for file in files:  # get the first psf file
        if "step4_pbcsetup.psf" in file:
            psf_file = os.path.join(active_dir_path, file)
            break
        elif "step3_pbcsetup.psf" in file:
            psf_file = os.path.join(active_dir_path, file)

    try:
        dcd = MDAnalysis.Universe(psf_file, dyna_file)
    except():
        print(f"!!!!!!!!!!!!!!!! DYNA FILE: {dyna_file}")
        print(f"!!!!!!!!!!!!!!!! PSF FILE: {psf_file}")
        print(f"!!!!!!!!!!!!!!!! PATH: {active_dir_path}")
        exit()

    iron_pos = dcd.trajectory[0]._pos[dcd.select_atoms("resid 1154 and type FE").ix[0]]
    return iron_pos


def get_N_position(active_dir_path):
    files = os.listdir(active_dir_path)
    psf_file = None
    dyna_file = None

    for file in files: # get the first dyna file
        if ".dcd" in file:
            dyna_file = os.path.join(active_dir_path, file)
            break
    for file in files: # get the first psf file
        if "step4_pbcsetup.psf" in file:
            psf_file = os.path.join(active_dir_path, file)
            break
        elif "step3_pbcsetup.psf" in file:
            psf_file = os.path.join(active_dir_path, file)

    if dyna_file is None or psf_file is None:
        print("No dyna or psf file found")
        print(active_dir_path)
        exit()

    dcd = MDAnalysis.Universe(psf_file, dyna_file)

    na_pos = dcd.trajectory[0]._pos[dcd.select_atoms("resid 1154 and name NA").ix[0]]
    nb_pos = dcd.trajectory[0]._pos[dcd.select_atoms("resid 1154 and name NB").ix[0]]
    nc_pos = dcd.trajectory[0]._pos[dcd.select_atoms("resid 1154 and name NC").ix[0]]
    nd_pos = dcd.trajectory[0]._pos[dcd.select_atoms("resid 1154 and name ND").ix[0]]

    pos_list = [na_pos, nb_pos, nc_pos, nd_pos]

    return pos_list


def calc_dist_to_plane(fe_pos, triplet):
    # Vectors in the plane
    AB = triplet[1] - triplet[0]
    AC = triplet[2] - triplet[0]

    # Normal vector
    normal = np.cross(AB, AC)

    # Plane equation: ax + by + cz + d = 0
    a, b, c = normal
    d = -np.dot(normal, triplet[0])

    #print(f"Plane equation: {a}x + {b}y + {c}z + {d} = 0")

    point_on_plane = triplet[0]
    point_fe_vect = fe_pos - point_on_plane
    distance = np.dot(normal, point_fe_vect) / np.linalg.norm(normal)

    # print(f"distance to plain: {distance}")
    return distance


def calc_avg_dist_to_plane(fe_pos, n_pos):
    # get all possible 3-point combinations for the plane of the four n atoms
    n_plane_triplets = [[n_pos[0], n_pos[1], n_pos[2]],
                        [n_pos[0], n_pos[1], n_pos[3]],
                        [n_pos[0], n_pos[2], n_pos[3]],
                        [n_pos[1], n_pos[2], n_pos[3]]]

    dist_list = []

    for triplet in n_plane_triplets:
        dist_list.append(calc_dist_to_plane(fe_pos, triplet))

    dist_list = np.array(dist_list)
    n_pos = np.array(n_pos)
    mean_pos = np.mean(n_pos, axis=0)
    # print(f"mean position of n atoms: {mean_pos}")
    # print(f"mean distance to plane: {np.linalg.norm(fe_pos-mean_pos)}")

    return np.mean(dist_list)


def psf_exists(_dir):
    files = os.listdir(_dir)
    for file in files:
        if "step4_pbcsetup.psf" in file or "step3_pbcsetup.psf":
            return True
    return False


def clearLog():
    if os.path.exists("distance_fe_plane_log.txt"):
        with open("distance_fe_plane_log.txt", "w") as f:
            f.write("")
        f.close()
        return


def log2file(text):
    if not os.path.exists("distance_fe_plane_log.txt"):
        with open("distance_fe_plane_log.txt", "w") as f:
            f.write(text)
        f.close()
        return

    with open("distance_fe_plane_log.txt", "a") as f:
        f.write(text)
    f.close()

def plot_temp_dir(curr_dir):
    sub_dirs = sorted(os.listdir(curr_dir))
    #print("processing directory " + curr_dir)

    possible_sub_dirs = np.arange(1, 51).tolist()

    for i in range(len(possible_sub_dirs)):
        possible_sub_dirs[i] = str(possible_sub_dirs[i])

    plane_distances = []

    for sub_dir in sub_dirs:
        if sub_dir in possible_sub_dirs:
            active_dir = os.path.join(curr_dir, sub_dir)
            files_in_sub_dir = os.listdir(active_dir)
            #print(f"processing subdir {sub_dir}")

            if "dyna0.dcd" in files_in_sub_dir and psf_exists(active_dir) or (not live and "dyna.dcd" in files_in_sub_dir):
                print(f"Processing directory {active_dir}")
                fe_pos = get_iron_position(active_dir)
                n_positions = get_N_position(active_dir)

                plane_dist = calc_avg_dist_to_plane(fe_pos, n_positions)
                plane_distances.append(plane_dist)
                #print(f"distance to plane: {plane_dist}")

                log2file(f"In directory {active_dir}: distance to plane: {plane_dist} \n")

    if len(plane_distances) > 0:
        #print(f"average distance to plane: {np.mean(plane_distances)}")
        # Plot histogram with improved appearance

        plt.hist(plane_distances, bins=15, edgecolor="black", alpha=0.75, color="royalblue", density=True)

        # Labels and title
        plt.xlabel("Distance Fe-Plane (Å)", fontsize=12)
        plt.ylabel("Probability", fontsize=12)
        plt.title("Histogram of Distances from N-plane to Fe", fontsize=14)

        # Grid for better readability
        plt.grid(axis="y", linestyle="--", alpha=0.7)
        if live:
            if "WITH_CO_V2" in curr_dir:
                plt.savefig("pictures/fe_plane_dist/CO/" + curr_dir.split("/")[-1] + ".png")
            else:
                plt.savefig("pictures/fe_plane_dist/NO_CO/" + curr_dir.split("/")[-1] + ".png")
        if not live:
            plt.show()


live = True
clearLog()

if live:
    cwd = os.getcwd()
else:
    cwd = os.getcwd() + "/../SimulationSetupAndRun"

sim_dirs = [cwd + "/WITH_CO_V2/", cwd + "/NO_CO_V2/"]
# temp_dirs = ["1K", "5K", "10K", "50K", "100K", "300K"]
temp_dirs = ["50K", "100K", "300K"]

if not live:
    temp_dirs = ["1K", "5K"]

total_dirs = []
for temp_dir in temp_dirs:
    total_dirs.append(sim_dirs[1] + temp_dir)
    total_dirs.append(sim_dirs[0] + temp_dir)

if not live:
    for _dir in total_dirs:
        plot_temp_dir(_dir)
else:
    with Pool(processes=len(total_dirs)) as pool:
        pool.map(plot_temp_dir, total_dirs)
