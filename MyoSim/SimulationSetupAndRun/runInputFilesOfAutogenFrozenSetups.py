import subprocess
import os

def get_subfolders():
    return os.listdir("../")

paths = get_subfolders()

path_prefixes = ["1K", "5K", "10K", "20K", "50K", "100K"]

folder_paths = []
start_read = 1
end_read = 2
WorWoCO = ["NO_CO", "WITH_CO"]

for co in WorWoCO:
    for tag in path_prefixes:
        for i in range(start_read, end_read):
            folder_paths.append(os.path.abspath(co + "/" + tag + "/" + str(i)))

print(folder_paths)

for path in folder_paths:
    result = subprocess.run(['sbatch', 'shell.sh'], capture_output=True, text=True, cwd=path)
    text = result.stdout

    with open("outfile_sh.out", "w") as file:
        file.write(text)
