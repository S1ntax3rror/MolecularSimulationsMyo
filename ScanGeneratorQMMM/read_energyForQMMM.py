import subprocess
import re
from glob import glob

import numpy as np
from matplotlib import pyplot as plt


output_file_path_prefix = "/home/kaeserj/PycharmProjects/CurveFitMorse/ScanGeneratorQMMM/outputMyoScan/output_orca_scan_"

input_file_list = glob(output_file_path_prefix + "*")
print(input_file_list)
print(len(input_file_list))
num_files = len(input_file_list)


energy_line_prefix = "FINAL SINGLE POINT ENERGY (QM/MM)"

energy_list = []

for i in range(num_files):
    output_file_path = output_file_path_prefix + f"{i:03d}"

    found = False

    with open(output_file_path, "r") as file:
        lines = file.readlines()
    file.close()

    for line in lines:
        if energy_line_prefix in line:
            line = float(line.split()[-1])
            energy_list.append((line))
            found = True

    if not found:
        energy_list.append(100000)



min_elem = min(energy_list)

for i in range(len(energy_list)):
    if energy_list[i] > 20000:
        energy_list[i] = 0
    else:
        energy_list[i] -= min_elem

setter = 0.7

while max(energy_list) > setter:
    energy_list[energy_list.index(max(energy_list))] = setter

print(max(energy_list))

distH2Fe = []
distH2 = []


with open("/home/kaeserj/PycharmProjects/CurveFitMorse/ScanGenerator/distancesH2AndFeH2ForQMMMAttempt1.txt", "r") as file:
    lines = file.readlines()

for line in lines:
    line = line.split()
    # print(line)
    distH2Fe.append(float(line[0]))
    distH2.append(float(line[1]))

print(len(distH2), len(distH2Fe), len(energy_list))

plt.tricontourf(distH2, distH2Fe, energy_list, 10, cmap="inferno_r")
plt.colorbar()
plt.ylabel("distance COM(H$_2$) to Fe")
plt.xlabel(r"bond distance H$_2$ ($\mathrm{\AA}$)")
plt.show()