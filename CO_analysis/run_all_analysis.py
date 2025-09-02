import os
import subprocess

main_dirs = ["NO_CO_V2", "NO_CO_V3", "WITH_CO_V2", "WITH_CO_V3"]

for main in main_dirs:
    for sub in os.listdir(main):
        sub_path = os.path.join(main, sub)
        if os.path.isdir(sub_path):
            try:
                subprocess.run(["sbatch", "run_analysis.sh"], cwd=sub_path)
                print(f"Submitted in {sub_path}")
            except subprocess.CalledProcessError as e:
                print(f"Failed in {sub_path}: {e}")
