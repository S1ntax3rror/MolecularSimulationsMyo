import os


def get_subdirs(path):
    files = os.listdir(path)
    dirs = []
    for file in files:
        if os.path.isdir(os.path.join(path, file)):
            dirs.append(os.path.join(path, file) + "\n")
    return dirs


def get_subdirs_to_use(path):
    files = os.listdir(path)
    dirs = []
    for file in files:
        if os.path.isdir(os.path.join(path, file)):
            dirs.append(os.path.join(path, file))
    return dirs


if __name__ == "__main__":
    cwd = os.getcwd()

    main_dirs = ["NO_CO_V2", "NO_CO_V3", "WITH_CO_V2", "WITH_CO_V3"]

    for dd in main_dirs:
        main_path = os.path.join(cwd, dd)
        temperature_dirs_ = get_subdirs_to_use(main_path)

        print(temperature_dirs_)
        for temperature_dir in temperature_dirs_:
            temp_dir_path = os.path.join(main_path, temperature_dir)
            print(temp_dir_path)
            dirs_ = get_subdirs(temp_dir_path)

            write_path = os.path.join(temp_dir_path, "subfolders.txt")
            print(f"writing {write_path}")
            with open(write_path, "w") as f:
                f.writelines(dirs_)