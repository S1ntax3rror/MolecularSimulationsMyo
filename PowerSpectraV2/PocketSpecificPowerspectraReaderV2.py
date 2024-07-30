import numpy as np

# Statistics
from statsmodels.tsa.stattools import acovf

# Matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

# Miscellaneous
import ase.units as units


def moving_average(data_set, periods=9):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, 'same')


def get_paths(prefix, suffix, num_h):
    paths = []
    for i in range(num_h):
        paths.append(prefix + str(i + 1) + suffix)
    return paths


#
#plot consts
#

# Morse MDCM
# x_lim_min = 4407.
# x_lim_max = 4527.

# CGenFF
x_lim_min = 4013.
x_lim_max = 4133.

normconst = 1.3

#
# replace distfile with the path to the H2 distances
#
# configure which data to include with the load_from
# options: MDCM, singleMDCM, noMDCM, MDCMmorse, singleMDCMmorse, noMDCMmorse
#


load_from = "noMDCM"
pocket_to_read = 5
all_pockets = True
use_consistent_size = True
frame_size = 20000  # only used if consistent size is activated
h2_to_show = [1, 2, 3, 4, 5]

load_from = load_from.split(" ")
pocket_libs = []

if "MDCM" in load_from:
    prefix = "MDCMWith5H2/h2distsV.H2-"
    suffix = ".dcd.1-10.npy"
    paths = get_paths(prefix, suffix, 5)
    pocket_libs.append("PocketLibrary/withMDCM_lib.txt")
if "noMDCM" in load_from:
    prefix = "NotMDCMWith5H2/h2distsV.H2-"
    suffix = ".dcd.1-90.npy"
    paths = get_paths(prefix, suffix, 5)
    pocket_libs.append("PocketLibrary/without_MDCM_lib.txt")
if "singleMDCM" in load_from:
    prefix = "SingleH2WithMDCM/h2distsV.SOLOH2-"
    suffix = ".dcd.1-10.npy"
    pocket_libs.append("PocketLibrary/soloMDCM_lib.txt")
    paths = get_paths(prefix, suffix, 2)
if "singleMDCMmorse" in load_from:
    prefix = "MDCMWithMorseSingleH2/h2distsV13.H2_"
    suffix = ".npy"
    pocket_libs.append("PocketLibrary/soloMDCM_Morse_lib.txt")
    paths = get_paths(prefix, suffix, 5)
if "MDCMmorse" in load_from:
    prefix = "MDCMmorse/h2distsV13.5_H2.H2-"
    suffix = ".dcd.1-38.npy"
    pocket_libs.append("PocketLibrary/MorseWithMDCM_lib.txt")
    paths = get_paths(prefix, suffix, 5)
if "noMDCMmorse" in load_from:
    prefix = "MDCMnoMorse/h2distsV12.5_H2.H2-"
    suffix = ".dcd.1-49.npy"
    pocket_libs.append("PocketLibrary/MorseNoMDCM_lib.txt")
    paths = get_paths(prefix, suffix, 5)

# contains lists containing num H2, pocket, start, end
pocket_info_matrix = []

for lib in pocket_libs:
    with open(lib, "r") as file:
        lines = file.readlines()

    for line in lines:
        line = line.split(" ")
        if len(line) == 4:
            if all_pockets:
                pocket_info_matrix.append([int(line[0]), int(line[1]), int(line[2]), int(line[3])])
            elif int(line[1]) == pocket_to_read:
                pocket_info_matrix.append([int(line[0]), int(line[1]), int(line[2]), int(line[3])])

freq = []
spec = []
ampl_h2 = []
pocket_number = []

cnt = -1


def helper(cnt, freq, spec, ampl_h2, pocket_number):
    cnt += 1
    pocket_number.append(data_list)
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
    h2_distances = h2_distances_loader[start_frame:end_frame]

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
    if (freq[cnt][1] - freq[cnt][0]) > favg:
        Nave = 5
    else:
        Nave = int(favg / (freq[cnt][1] - freq[cnt][0]))
    print("moving avg searchradius):", str(Nave))
    spec[cnt] = moving_average(spec[cnt], Nave)

    # Get CN vibration amplitude
    ampl_h2.append(np.max(spec[cnt][np.logical_and(freq[cnt] > 3500., freq[cnt] < 4500.)]))
    return cnt, freq, spec, ampl_h2, pocket_number


def helperForCompleteSpectra(file):
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

    h2_distances = np.load(file)

    num_timesteps = len(h2_distances)
    time = np.arange(0, num_timesteps) * timedt

    # -----------------------------
    # Predict IR spectra
    # -----------------------------

    # Frequency array
    nf = len(time)
    dt = time[1] - time[0]
    Nfrq = int(nf / 2) + 1
    freq2 = np.arange(Nfrq) / float(nf) / dt * jiffy

    # Compute IR spectra
    acv = acovf(h2_distances, fft=True)

    acv = acv * np.blackman(nf)
    spec2 = np.abs(np.fft.rfftn(acv))

    beta = 1.0 / (kB_Ha_K * au2cminv) / T
    spec2 = spec2 * freq2 * (1 - np.exp(-freq2 * beta))

    # -----------------------------
    # Plot IR spectra
    # -----------------------------

    # Apply moving average

    favg = 5.0
    if (freq2[1] - freq2[0]) > favg:
        Nave = 5
    else:
        Nave = min(int(favg / (freq2[1] - freq2[0])), 50)
    print("moving avg searchradius):", str(Nave))
    spec2 = moving_average(spec2, Nave)

    # Get CN vibration amplitude
    ampl_h2 = spec2[np.logical_and(freq2 > 3500., freq2 < 4500.)]
    return spec2, freq2, ampl_h2


for data_list in pocket_info_matrix:
    distfile = paths[data_list[0] - 1]

    start_frame = data_list[2]
    end_frame = data_list[3]

    if use_consistent_size:
        if end_frame - start_frame < frame_size:
            continue
        start_frame = data_list[2]
        end_frame = start_frame + frame_size

        while end_frame <= data_list[3]:
            cnt, freq, spec, ampl_h2, pocket_number = helper(cnt, freq, spec, ampl_h2, pocket_number)

            start_frame = end_frame
            end_frame += frame_size

    else:
        cnt, freq, spec, ampl_h2, pocket_number = helper(cnt, freq, spec, ampl_h2, pocket_number)

spec_avg = np.zeros_like(spec[0])

total_spec = []
for i in range(len(h2_to_show)):
    one_spec, freq2, ampl_h2 = helperForCompleteSpectra(paths[i])
    if len(total_spec) != 0:
        total_spec += np.array(one_spec)
    else:
        total_spec, freq2, ampl_h2 = helperForCompleteSpectra(paths[i])
        total_spec = np.array(total_spec)
# print(total_spec.shape)
# print(len(one_spec))
# print(len(freq2))
# np.save(spec_avg, "spec_avg" + load_from + "all_pockets")

# for i in range(len(freq)):
#     if pocket_number[i][0] in h2_to_show:
#         spec_avg += spec[i]
#         numbers = [freq[i], spec[i] / ampl_h2[i]]
#         #for
#         #plt.plot(*numbers, label=str(pocket_number[i]), linewidth=2.0)


SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

fig = plt.figure(figsize=(8, 8))
axs = fig.add_subplot(1, 1, 1)

cnt = 0

color_scheme = [
    'b', 'r', 'g', 'purple', 'orange', 'magenta', 'brown', 'darkblue',
    'darkred', 'darkgreen', 'darkgrey', 'olive']

pockets_to_display = [2, 4, 5, 7]

spec_avg_list = []

for pocket_i in range(9):
    spec_avg = np.zeros_like(spec[0])
    available = False
    for i in range(len(freq)):
        if pocket_number[i][1] == pocket_i + 1:
            spec_avg += spec[i]
            available = True
    if available and pocket_i + 1 in pockets_to_display:
        cnt += 1
        spec_avg_list.append(spec_avg)
        ampl_avg = np.max(spec_avg[freq[0] > 3500.])
        plt.plot(freq[0], spec_avg / (ampl_avg * normconst) + 5 / 1.8 - cnt / 1.8, label=str(pocket_i + 1),
                 linewidth=2.0, color=color_scheme[pocket_i])
#        print(pocket_number[i][0], pocket_number[i][1], pocket_number[i][2], pocket_number[i][3])

np_spec_avg = np.array(spec_avg_list)
np.savez("spec_avg_save", spec_avg=np_spec_avg, freq=freq[0])

ampl_avg = np.max(total_spec[freq2 > 3500.])
plt.plot(freq2, total_spec / (ampl_avg * normconst), 'k', label=str("total"), linewidth=4.0)

tbox = TextArea(
    'A',
    textprops=dict(
        color='k', fontsize=35, ha='center', va='center')
)
anchored_tbox = AnchoredOffsetbox(
    loc="upper right", child=tbox, pad=0., frameon=False,
    bbox_to_anchor=(0.97, 0.97),
    bbox_transform=axs.transAxes, borderpad=0.)
axs.add_artist(anchored_tbox)

plt.xlim(x_lim_min, x_lim_max)

# plt.title("Powerspectra Average", fontsize=25)

plt.yticks([])

if load_from[0] == "MDCMmorse":
    loader = np.load("middle_mdcm_morse.npz")
    pos = loader["pos"]
    cnt = 0
    for tup in pos:
        pocket_i = pockets_to_display[cnt]
        cnt += 1
        print(tup)
        plt.scatter(tup[0], 5 / 1.8 - (cnt-1) / 1.8, color=color_scheme[pocket_i], marker="o")
        # plt.scatter(tup[1], 5 / 1.8 - (cnt-1) / 1.8,color=color_scheme[pocket_i], marker="o", alpha=0.6)

if load_from[0] == "noMDCM":
    loader = np.load("middle_NoMDCM.npz")
    pos = loader["pos"]
    cnt = 0
    for tup in pos:
        pocket_i = pockets_to_display[cnt]
        cnt += 1
        print(tup)
        plt.scatter(tup[0], 5 / 1.8 - (cnt-1) / 1.8,color=color_scheme[pocket_i], marker="o")
        # plt.scatter(tup[1], 5 / 1.8 - (cnt-1) / 1.8,color=color_scheme[pocket_i], marker="o", alpha=0.6)



plt.xlabel("Frequency (cm$^{-1}$)")
plt.ylabel("Intensity (arb. units)")

plt.legend(loc='upper left', title='H$_2$ in Pocket', framealpha=1.0)

plt.savefig(dpi=200, fname=str(load_from[0]), )

plt.show()
