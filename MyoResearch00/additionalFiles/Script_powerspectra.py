# Basics
import os
import sys

import MDAnalysis
import numpy as np
from glob import glob

# Statistics
from statsmodels.tsa.stattools import acovf

# Matplotlib
import matplotlib
import matplotlib.pyplot as plt

# Miscellaneous
import ase.units as units


def read_pocket_list(res_list):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd.select_atoms(res).ix[0])
    return np.array(pocket_indices)


# -----------------------------
# Parameters
# -----------------------------

# Temperatures (K)
T = 300.0

# Time for speed of light in vacuum to travel 1 cm in ps
jiffy = 0.01 / units._c * 1e12
# 33.3564 [cm/ps]

# Boltzman constant in Hartree/Kelvin
kB_Ha_K = units.kB / units.Hartree
# 3.16681e-06 [Hartree/K]

# Conversion factor from Hartree (au) to inverse centimeter cm**-1
# 1 Hartree = [Hartree to Joule) / (planck constant [Js] * speed of light [m/s]) * (1/m to 1/cm)] cm**-1
au2cminv = 1.0 * units.Hartree / units.J / (units._hplanck * units._c) * 1.e-2
# 219474.63 [1/Hartree/cm]

# -----------------------------
# Get data
# -----------------------------

split_dcd = ['.', 1]  # split a each '.' and take element [1] in filename_dcd

path_psf = '/Data/dcdForH2Bondwiggle/step3_pbcsetup.psf'

datadir = '/Data/dcdForH2Bondwiggle'

file_dcd = '/home/kaeserj/PycharmProjects/CurveFitMorse/Data/dcdForH2Bondwiggle/dyna1.dcd'

num_timesteps = 0
timedt = 0.00025  # Time / frame in ps

# ... Get H2 bond distances in Angstrom
H2_atoms = ["resname H2 and resid 156 and name H1", "resname H2 and resid 156 and name H2"]
h2_distances = []

# Open dcd file
dcd = MDAnalysis.Universe(path_psf, file_dcd)

# Get atom indices for H2
H2_arr = read_pocket_list(H2_atoms)

Nframes = len(dcd.trajectory)
Nskip = int(dcd.trajectory.skip_timestep)
dt = np.round(
    float(dcd.trajectory._ts_kwargs['dt']), decimals=8) / Nskip

atoms = np.array([ai for ai in dcd._topology.names.values])
Natoms = len(atoms)

distfile = 'h2dists.npy'

if not os.path.exists(distfile):

    for jj, frame in enumerate(dcd.trajectory):
        num_timesteps += 1
        # Get updated H2 position
        positions_H2 = frame._pos[H2_arr]
        h2_dist = np.linalg.norm((positions_H2[0] - positions_H2[1]))
        h2_distances.append(h2_dist)

    # ... Get time sequence [0.00025, 0.0005, 0.00075, ...] in ps
    time = np.arange(0, num_timesteps) * timedt
    h2_distances = np.array(h2_distances)

    np.save(distfile, h2_distances)

else:

    h2_distances = np.load(distfile)
    num_timesteps = len(h2_distances)
    time = np.arange(0, num_timesteps) * timedt

# -----------------------------
# Predict IR spectra
# -----------------------------

# Frequency array
nf = len(time)
dt = time[1] - time[0]
Nfrq = int(nf / 2) + 1
freq = np.arange(Nfrq) / float(nf) / dt * jiffy

# Compute IR spectra
acv = acovf(h2_distances, fft=True)

acv = acv * np.blackman(nf)
spec = np.abs(np.fft.rfftn(acv))

beta = 1.0 / (kB_Ha_K * au2cminv) / T
spec = spec * freq * (1 - np.exp(-freq * beta))


# -----------------------------
# Plot IR spectra
# -----------------------------

# Apply moving average
def moving_average(data_set, periods=9):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, 'same')


favg = 1.0
Nave = int(favg / (freq[1] - freq[0]))
# spec = moving_average(spec, Nave)

# Get CN vibration amplitude
ampl_h2 = np.max(spec[np.logical_and(freq > 3500., freq < 4500.)])

plt.plot(freq, spec / ampl_h2)

plt.xlim(500., 4500)

plt.show()
