import subprocess
import os


def get_subfolders():
    return os.listdir("../")


paths = get_subfolders()

path_prefixes = ["1K", "5K", "10K", "20K", "50K", "100K"]

dir_list = []

for path in paths:
    form_path = path.split("_")[0]
    if any(p == form_path for p in path_prefixes):
        dir_list.append(path)

sh_file_paths = []
for path in dir_list:
    sh_file_paths.append(os.path.abspath(path))
print(sh_file_paths)

for path in sh_file_paths:
    result = subprocess.run(['sbatch', 'shell.sh'], capture_output=True, text=True, cwd=path)
    text = result.stdout

    with open("outfile_sh.out", "w") as file:
        file.write(text)
