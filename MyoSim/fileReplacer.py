import glob
import os


def get_subfolders():
    return os.listdir("./")


def replace_line_in_file(path_, filename, line):
    act_path = "./" + path_[0] + "/" + filename
    print(act_path)

    with open(act_path, "r") as prod:
        p_lines = prod.readlines()
        for i, line_ in enumerate(p_lines):
            print(i, line_)
        p_lines[line] = "set temp = " + path_[1] + "\n"
        for i, line_ in enumerate(p_lines):
            print(i, line_)
    prod.close()
    with open(act_path, "w") as file:
        file.writelines(p_lines)
    file.close()


paths = get_subfolders()

path_prefixes = ["1K", "5K", "10K", "20K", "50K", "100K"]

dir_list = []

for path in paths:
    form_path = path.split("_")[0]
    if any(p == form_path for p in path_prefixes):
        dir_list.append((path, form_path.strip("K")))


for path in dir_list:
    replace_line_in_file(path, "prod22.inp", 88)
    replace_line_in_file(path, "step4_equilibration.inp", 127)



    # with open(path+"step4_equilibration.inp", "w") as step4:
    #     s_lines = step4.readlines()
    # step4.close()
