import os

import MDAnalysis
import numpy as np
from matplotlib import pyplot as plt


def get_coordinates(path_psf_, dyna_files_, atom="O"):
    """
    :param path_psf_: path to PSF file
    :param dyna_files_: list of MD dynamics files
    :param atom: Selection parameter [C or O]
    :return: list of coordinates, dcd refference
    """
    if atom != "O" and atom != "C":
        atom="O"

    dcd = MDAnalysis.Universe(path_psf_, dyna_files_)
    coords = np.zeros((len(dcd.trajectory), 3))
    CO = dcd.select_atoms(f"resname CO and name {atom}")

    for i, ts in enumerate(dcd.trajectory):
        coords[i] = CO.positions

    return coords, dcd


def get_dcds_in_order(folder_path):
    """
    Returns all dcd files in a folder sorted correcly as long as they follow the naming: dyna{i}.dcd
    :param folder_path:
    :return:
    """
    files = os.listdir(folder_path)
    dcd_files_ = []
    for i in range(1000):
        if f"dyna{i}.dcd" in files:
            dcd_files_.append(f"dyna{i}.dcd")
    return dcd_files_


def get_distances_list_FE_coords(FE_position, Coordinates_list):
    """
    Takes in the iron position and the coordinates list. Returns the distance between the iron and the coordinates
    :param FE_position:
    :param Coordinates_list:
    :return:
    """
    return np.linalg.norm(Coordinates_list - FE_position, axis=1)


def get_iron_position(dcd_universe):
    """
    Returns the iron position given a valid MDUniverse
    :param dcd_universe:
    :return:
    """
    FE = dcd_universe.select_atoms("resid 1154 and type FE")
    iron_pos = [0, 0, 0]
    for i, ts in enumerate(dcd_universe.trajectory):
        iron_pos = FE.positions
        break
    return iron_pos


def plot_distances(distances_, bins=30):
    """
    Takes a list of distances and plots the probability distibution of the distances
    :param distances_:
    :param bins:
    :return:
    """
    plt.hist(distances_, bins=bins, edgecolor="black", alpha=0.75, color="royalblue", density=True)

    # Labels and title
    plt.xlabel("Distance to Fe (Å)", fontsize=12)
    plt.ylabel("Probability", fontsize=12)
    plt.title("Histogram of Distances to Fe", fontsize=14)

    # Grid for better readability
    plt.grid(axis="y", linestyle="--", alpha=0.7)


    cwd_ = os.getcwd()
    if "/home/jiri/PycharmProjects/" in cwd_:
        plt.savefig("pic.png")
        plt.show()
    else:
        cwd_ = cwd_.split("/")
        plt.savefig(f"FE_CO_distance_plot{cwd_[-2]}_{cwd_[-1]}.png")


if __name__ == "__main__":
    psf_files = "step4_pbcsetup.psf"
    cwd = os.getcwd()
    if os.path.exists("distances.npy"):
        distances = np.load("distances.npy")
    else:
        dcd_files = get_dcds_in_order(cwd)
        coords, dcd_univ = get_coordinates(psf_files, dcd_files)
        fe_coordinates = get_iron_position(dcd_univ)
        distances = get_distances_list_FE_coords(fe_coordinates, coords)
        np.save("distances.npy", distances)
    plot_distances(distances)
