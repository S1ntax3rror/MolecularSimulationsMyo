import numpy as np
import matplotlib.pyplot as plt
from ase import io
from scipy.optimize import curve_fit

f = open("dimorse_test.dat", "r")
numOfLines = 100

energy = [0] * numOfLines
dist1 = [0] * numOfLines
dist2 = [0] * numOfLines

for i in range(numOfLines):
    dist1[i], dist2[i], energy[i] = f.readline().split(" ")

for i in range(len(dist2)):
    dist2[i] = ''.join(dist2[i].split())[:-1].upper()

    dist2[i] = float(dist2[i])
    dist1[i] = float(dist1[i])
    energy[i] = float(energy[i])

dist1 = np.array(dist1)
dist2 = np.array(dist2)


def morse(d, re, De, beta):
    d1 = d[0, :]
    d2 = d[1, :]
    return De * (1. - np.exp(-beta * (d1 - re))) ** 2 + De * (1. - np.exp(-beta * (d2 - re))) ** 2


xdata = np.concatenate((dist1.reshape(1, -1), dist2.reshape(1, -1)), axis=0)
popt, pcov = curve_fit(morse, xdata, energy, maxfev=40000000, p0=[7.81145664, 1.5, 1])

yfit = morse(xdata, popt[0], popt[1], popt[2])

print(popt)

plt.plot(xdata[0], energy, "ro")
plt.plot(xdata[0], yfit)

plt.show()
