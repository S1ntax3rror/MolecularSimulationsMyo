import numpy as np
import matplotlib.pyplot as plt
from ase import io
from scipy.optimize import curve_fit

f = open("energyData2DScan.dat", "r")
numOfLines = 494

energy = np.zeros(numOfLines, dtype=float)
distH2Fe = np.zeros(numOfLines, dtype=float)
distH2 = np.zeros(numOfLines, dtype=float)

"""
with open("energyData.dat", "r") as f:
    lines = f.readlines()
for i, line in enumerate(lines):
    if not len(line):
        break
    height[i], energy[i] = line.split(" ")
"""

for i in range(numOfLines):
    # print(f.readline().split(" "))
    distH2Fe[i], distH2[i], energy[i] = np.array(f.readline().split(), dtype=float)

f.close()

energy = energy - np.min(energy)

#for i in range(0, numOfLines):
#    print(str(distH2Fe[i]).strip("\n") + "    " + str(distH2[i]).strip("\n") + "    " + str(energy[i]).strip("\n"))

plt.tricontourf(distH2*2, distH2Fe, energy, 10, cmap="inferno_r")
plt.colorbar()
plt.ylabel("distance COM(H$_2$) to Fe")
plt.xlabel(r"bond distance H$_2$ ($\mathrm{\AA}$)")
plt.show()