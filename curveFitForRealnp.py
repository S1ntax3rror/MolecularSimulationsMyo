import numpy as np
import matplotlib.pyplot as plt
from ase import io
from scipy.optimize import curve_fit


f = open("energyData.dat", "r")
numOfLines = 21

energy = [0] * numOfLines
height = [0] * numOfLines

dist1 = [0] * numOfLines
dist2 = [0] * numOfLines
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
    height[i], energy[i] = f.readline().split(" ")

f.close()
energy = np.array(energy, dtype=float)

systems = io.read("heme_scan_fe_h2_v2.xyz")
for i in range(numOfLines):
    systems = io.read("heme_scan_fe_h2_v2.xyz", index=i)
    dist1[i] = systems.get_distance(17, 77)
    dist2[i] = systems.get_distance(17, 78)

dist1 = np.array(dist1)
dist2 = np.array(dist2)


def morse(d, De, re, beta, y):
    d1 = d[0, :]
    d2 = d[1, :]
    return De * (1. - np.exp(-beta * (d1 - re))) ** 2 + De * (1. - np.exp(-beta * (d2 - re))) ** 2 + y


xdata = np.concatenate((dist1.reshape(1, -1), dist2.reshape(1, -1)), axis=0)

popt, pcov = curve_fit(morse, xdata, energy, maxfev=40000000, p0=[7, 1.65, 1, -0.215])

yfit = morse(xdata, popt[0], popt[1], popt[2], popt[3])

print(popt)

plt.plot(xdata[0], energy, "ro")
plt.plot(xdata[0], yfit)

plt.show()
