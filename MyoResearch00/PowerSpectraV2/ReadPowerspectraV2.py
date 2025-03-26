import numpy as np

# Statistics
from statsmodels.tsa.stattools import acovf

# Matplotlib
import matplotlib.pyplot as plt

# Miscellaneous
import ase.units as units


def moving_average(data_set, periods=9):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, 'same')


#
# replace distfile with the path to the H2 distances
#


distfile = "MDCMWith5H2/h2distsV.H2-2.dcd.1-10.npy"

start_frame1 = 0
end_frame1 = 20000

start_frame2 = 00000
end_frame2 = 250000



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
# 1 Hartree = [Hartree to Joule] / (planck constant [Js] * speed of light [m/s]) * (1/m to 1/cm)] cm**-1
au2cminv = 1.0 * units.Hartree / units.J / (units._hplanck * units._c) * 1.e-2
# 219474.63 [1/Hartree/cm]


timedt = 0.001  # Time / frame in ps


h2_distances_loader = np.load(distfile)
h2_distances1 = h2_distances_loader[start_frame1:end_frame1]
h2_distances2 = h2_distances_loader[start_frame2:end_frame2]
h2_dist_arr = [h2_distances1, h2_distances2]


freq = []
spec = []
ampl_h2 = []

cnt = -1
for h2_distances in h2_dist_arr:
    cnt += 1

    num_timesteps = len(h2_distances)
    time = np.arange(0, num_timesteps) * timedt

    # -----------------------------
    # Predict IR spectra
    # -----------------------------

    # Frequency array
    nf = len(time)
    dt = time[1] - time[0]
    Nfrq = int(nf / 2) + 1
    freq.append(np.arange(Nfrq) / float(nf) / dt * jiffy)

    # Compute IR spectra
    acv = acovf(h2_distances, fft=True)

    acv = acv * np.blackman(nf)
    spec.append(np.abs(np.fft.rfftn(acv)))

    beta = 1.0 / (kB_Ha_K * au2cminv) / T
    spec[cnt] = spec[cnt] * freq[cnt] * (1 - np.exp(-freq[cnt] * beta))

    # -----------------------------
    # Plot IR spectra
    # -----------------------------

    # Apply moving average

    favg = 2.0
    Nave = int(favg / (freq[cnt][1] - freq[cnt][0]))
    print("moving avg searchradius):", str(Nave))
    spec[cnt] = moving_average(spec[cnt], Nave)

    # Get CN vibration amplitude
    ampl_h2.append(np.max(spec[cnt][np.logical_and(freq[cnt] > 3500., freq[cnt] < 4500.)]))


for i in range(len(freq)):
    plt.plot(freq[i], spec[i] / ampl_h2[i], label="pocket " + str(i))


plt.xlim(3000., 4500)

plt.title("h2specV11.H2-" + "1")

plt.legend(loc='upper left')

plt.show()

