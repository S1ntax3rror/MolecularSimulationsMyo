import os
from distutils.dir_util import copy_tree
from os import mkdir
import shutil
from multiprocessing import Pool

import numpy as np

from getPocketCenterOfPDB import add_co_at_center_of_pocket_5, get_coord_of_xe5, get_coord_of_CO
from getCenterBetweenXE5_4 import get_xe5_xe4_center


def replace_line_in_file(path_, line, replace_text):
    with open(path_, "r") as file:
        p_lines = file.readlines()
        p_lines[line] = replace_text + "\n"
    file.close()

    with open(path_, "w") as file:
        file.writelines(p_lines)
    file.close()


def generate_temperature_folder_without_co(tag_):

    basepath = no_co_sim_folder_path + "/" + tag_ + "/"
    all_folders_at_tag_ = os.listdir(basepath.rstrip("/"))

    for dyna_file_suffix_ in range(start_read, end_read):
        if str(dyna_file_suffix_) in all_folders_at_tag_:
            print("Skipping generation of folder", dyna_file_suffix_)
            continue

        print("generating folder without co", tag_, dyna_file_suffix_)

        i = tags.index(tag_)
        print(i)

        mkdir(basepath + str(dyna_file_suffix_))
        copy_tree(setup_no_co_path, basepath + str(dyna_file_suffix_))

        path_dyna_crd = r_folders_wo_co[i] + "/dyna" + str(dyna_file_suffix_) + ".crd"
        path_dyna_res = r_folders_wo_co[i] + "/dyna" + str(dyna_file_suffix_) + ".res"
        path_dyna_pdb = r_folders_wo_co[i] + "/dyna" + str(dyna_file_suffix_) + ".pdb"
        path_dyna_psf = r_folders_wo_co[i] + "/step3_pbcsetup.psf"
        shutil.copyfile(path_dyna_res, basepath + str(dyna_file_suffix_) + "/dyna.res")
        shutil.copyfile(path_dyna_crd, basepath + str(dyna_file_suffix_) + "/dyna.crd")
        if live:
            shutil.copyfile(path_dyna_psf, basepath + str(dyna_file_suffix_) + "/step3_pbcsetup.psf")

        center_XE5_XE4 = get_xe5_xe4_center(path_dyna_pdb)
        add_co_at_center_of_pocket_5(path_dyna_pdb, basepath + str(dyna_file_suffix_) + "/co.crd")

        path_step_5_inp = basepath + str(dyna_file_suffix_) + "/prod22.inp"
        path_step_4_inp = basepath + str(dyna_file_suffix_) + "/step4_equilibration.inp"
        print(center_XE5_XE4)
        replace_line_in_file(path_step_4_inp, 82,
                             "XREF " + str(center_XE5_XE4[0]) + " YREF " + str(center_XE5_XE4[1]) + " ZREF " + str(
                                 center_XE5_XE4[2]) + " -")
        replace_line_in_file(path_step_5_inp, 82,
                             "XREF " + str(center_XE5_XE4[0]) + " YREF " + str(center_XE5_XE4[1]) + " ZREF " + str(
                                 center_XE5_XE4[2]) + " -")

        if i < 6:
            replace_line_in_file(path_step_5_inp, 113, "set temp = " + str(tag.replace("K", "")))
            replace_line_in_file(path_step_4_inp, 98, "set temp = " + str(tag.replace("K", "")))

        else:
            replace_line_in_file(path_step_5_inp, 113, "set temp = " + str("303.15"))
            replace_line_in_file(path_step_4_inp, 98, "set temp = " + str("303.15"))


def generate_temperature_folder_with_co(tag_):
    basepath = with_co_sim_folder_path + "/" + tag_ + "/"
    all_folders_at_tag = os.listdir(basepath.rstrip("/"))

    i = tags.index(tag_)

    for dyna_file_suffix in range(start_read, end_read):
        if str(dyna_file_suffix) in all_folders_at_tag:
            print("Skipping generation of folder", dyna_file_suffix)
            continue

        mkdir(basepath + str(dyna_file_suffix))
        copy_tree(setup_co_path, basepath + str(dyna_file_suffix))

        print("generating folder with co", tag_, dyna_file_suffix)
        path_dyna_crd = r_folders_with_co[i] + "/dyna" + str(dyna_file_suffix) + ".crd"
        path_dyna_res = r_folders_with_co[i] + "/dyna" + str(dyna_file_suffix) + ".res"
        path_dyna_pdb = r_folders_with_co[i] + "/dyna" + str(dyna_file_suffix) + ".pdb"
        path_dyna_psf = r_folders_with_co[i] + "/step3_pbcsetup.psf"

        shutil.copyfile(path_dyna_res, basepath + str(dyna_file_suffix) + "/dyna.res")
        shutil.copyfile(path_dyna_crd, basepath + str(dyna_file_suffix) + "/dyna.crd")
        if live:
            shutil.copyfile(path_dyna_psf, basepath + str(dyna_file_suffix) + "/step3_pbcsetup.psf")

        path_step_5_inp = basepath + str(dyna_file_suffix) + "/prod22.inp"
        replace_line_in_file(path_step_5_inp, 113, "set temp = " + str(tag.replace("K", "")))

        pocket_center_coords = get_coord_of_xe5(path_dyna_pdb)

        CO_coords = get_coord_of_CO(path_dyna_pdb)
        dist_center_CO = np.linalg.norm(pocket_center_coords - CO_coords)
        print(dist_center_CO)

        if dist_center_CO > 5:
            print("CO escaped for " + tag + "dyna" + str(dyna_file_suffix))

# used to make switching between live and test system more easy
live = True

start_read = 1
end_read = 51

# init paths
setup_co_path = "SetupFiles_Frozen_init_with_co"
setup_no_co_path = "SetupFiles_Frozen_init_without_co"

frozen_myo_generator_path = "/cluster/data/toepfer/Project_MyoSim/GenerateFrozenMyo"
no_co_sim_folder_path = "NO_CO_V3"
with_co_sim_folder_path = "WITH_CO_V3"

if not live:
    start_read = 1
    end_read = 3
    frozen_myo_generator_path = "frozen_myo_strucutre_data"

read_from_folders = os.listdir(frozen_myo_generator_path)
filter = ["getHeatedMyo_StartedWithCO", "getHeatedMyo_V2_with_FE-HETA"]

if not live:
    filter = ["withCO", "woutCO"]

r_folders_with_co = []
r_folders_wo_co = []

for item in read_from_folders:
    if filter[0] in item:
        r_folders_with_co.append(frozen_myo_generator_path + "/" + item)
    elif filter[1] in item:
        r_folders_wo_co.append(frozen_myo_generator_path + "/" + item)
r_folders_wo_co = sorted(r_folders_wo_co)
r_folders_with_co = sorted(r_folders_with_co)
print(r_folders_with_co)
print(r_folders_wo_co)

tags = ["1K", "5K", "10K", "20K", "50K", "100K", "300K"]
tags = ["1K", "5K", "10K", "20K"]
if not live:
    tags = ["100K", "300K"]

write_to_folders_w_co = os.listdir(with_co_sim_folder_path)
write_to_folders_wo_co = os.listdir(no_co_sim_folder_path)

# if not live:
#     for folder in write_to_folders_wo_co:
#         os.removedirs(with_co_sim_folder_path + "/" + folder)
#     for folder in write_to_folders_w_co:
#         os.removedirs(no_co_sim_folder_path + "/" + folder)

# create sim folders
if len(write_to_folders_w_co) == 0:
    for tag in tags:
        os.mkdir(with_co_sim_folder_path + "/" + tag)
else:
    if len(write_to_folders_w_co) != 7:
        print("ERROR not enough folders defined")
        exit()
if len(write_to_folders_wo_co) == 0:
    for tag in tags:
        os.mkdir(no_co_sim_folder_path + "/" + tag)
else:
    if len(write_to_folders_wo_co) != 7:
        print("ERROR not enough folders defined")
        exit()

# # gen setups for sim that already had CO in it
# for tag in tags:
#     generate_temperature_folder_with_co(tag)
#
# # gen setups for sim that don't have CO in them yet
# for tag in tags:
#     generate_temperature_folder_without_co(tag)

if __name__ == "__main__":
    with Pool(processes=len(tags)) as pool:
        pool.map(generate_temperature_folder_with_co, tags)

    with Pool(processes=len(tags)) as pool:
        pool.map(generate_temperature_folder_without_co, tags)


