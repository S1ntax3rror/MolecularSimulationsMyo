import numpy as np
import matplotlib.pyplot as plt
from ase import io
from scipy.optimize import curve_fit

f = open("morseH2.dat", "r")
numOfLines = 100

energy = [0] * numOfLines
dist1 = [0] * numOfLines
for i in range(numOfLines):
    dist1[i], energy[i] = f.readline().split(" ")

for i in range(len(dist1)):

    dist1[i] = float(dist1[i])
    energy[i] = float(energy[i])

dist1 = np.array(dist1)
energy = np.array(energy)
energy -= min(energy)

def morse(d, De, re, beta, shift):
    d1 = d
    return De * (1. - np.exp(-beta * (d1 - re))) ** 2 + shift


xdata = dist1
popt, pcov = curve_fit(morse, xdata, energy, maxfev=40000000, p0=[5, 0.7, 0.7, -1.0])

yfit = morse(xdata, popt[0], popt[1], popt[2], popt[3])

print(popt)

plt.plot(xdata, energy, "ro")
plt.plot(xdata, yfit)

plt.show()

