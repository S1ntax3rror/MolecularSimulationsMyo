import os
import shutil
from os import mkdir
from distutils.dir_util import copy_tree


def generate_all_permutations_unbound(active_directory, base_setup_directory):
    line1 = "patch ULIG COA    1 HETA 142 PROA 87 noangle nodihedral"
    line2 = "patch ULIG COB    1 HETB 147 PROB 92 noangle nodihedral"
    line3 = "patch ULIG COC    1 HETC 142 PROC 87 noangle nodihedral"
    line4 = "patch ULIG COD    1 HETD 147 PROD 92 noangle nodihedral"

    for bond_1 in range(2):
        for bond_2 in range(2):
            for bond_3 in range(2):
                for bond_4 in range(2):
                    print("setting up folder for permuation:: ", bond_1, bond_2, bond_3, bond_4)
                    current_dir = active_directory + "/" + "permutation" + str(bond_1) + str(bond_2) + str(
                        bond_3) + str(bond_4)
                    mkdir(current_dir)
                    copy_tree(base_setup_directory, current_dir)
                    path_step6 = current_dir + "/" + "gpu_step6_manipulation.inp"
                    with open(path_step6, "r") as f:
                        text = f.readlines()
                        text_copy = text.copy()
                        for i, line in enumerate(text):
                            if line1 in line and bond_1 == 1:
                                text_copy[i] = "!" + line1 + "\n"
                            elif line2 in line and bond_2 == 1:
                                text_copy[i] = "!" + line2 + "\n"
                            elif line3 in line and bond_3 == 1:
                                text_copy[i] = "!" + line3 + "\n"
                            elif line4 in line and bond_4 == 1:
                                text_copy[i] = "!" + line4 + "\n"
                        f.close()
                    with open(path_step6, "w") as f:
                        f.writelines(text_copy)
                        f.close()


if __name__ == "__main__":

    overwrite = False
    live = True

    cwd = os.getcwd()

    if "cluster" in cwd:
        live = True

    hemo_base_path = cwd + "/" + "HemoBase"
    path_setup_destination = "all_permutations_unbound"

    if live:
        hemo_base_path = cwd + "/" + "Kai_setup"

    if os.path.exists(cwd + "/" + path_setup_destination) and (not overwrite and not live):
        print("EITHER SET OVERWRITE TRUE OR CHANGE PATH")
        exit()
    elif os.path.exists(cwd + "/" + path_setup_destination):
        shutil.rmtree(cwd + "/" + path_setup_destination)

    mkdir(cwd + "/" + path_setup_destination)
    active_dir = cwd + "/" + path_setup_destination

    generate_all_permutations_unbound(active_dir, hemo_base_path)
