import subprocess
import os


def replace_line_in_file(file, old_line, new_line):
    with open(file, "r") as f:
        lines = f.readlines()
        lines_copy = lines.copy()
        for i, line in enumerate(lines):
            if old_line in line:
                lines_copy[i] = new_line
        f.close()
    with open(file, "w") as f:
        f.writelines(lines_copy)
        f.close()


if __name__ == "__main__":
    cwd = os.getcwd()
    active_dir = cwd + "/" + "all_permutations_unbound/"

    paths = os.listdir(active_dir)

    sh_file_paths = []

    for path_ in paths:
        run_path = os.path.join(active_dir, path_) + "/"
        sh_file_paths.append(run_path)

    print(sh_file_paths)
    good_gpus = ["gpu04", "gpu25", "gpu06", "gpu26", "gpu03", "gpu12", "gpu13", "gpu18", "gpu01", "gpu02"]

    for path in sh_file_paths:
        result = subprocess.run(['sbatch', 'run.sh'], capture_output=True, text=True, cwd=path)
        text = result.stdout
        text_err = result.stderr

        with open("outfile_sh.out", "w") as file:
            file.write(text)
            file.write(text_err)
