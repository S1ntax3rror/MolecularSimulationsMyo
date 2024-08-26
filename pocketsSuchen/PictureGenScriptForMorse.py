import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
from statistics import mode

prefix = "pocket_arrays/"
suffix = 'pocket_arrays_v14.5_H2.DCD_files.1-37.H2_num2.npz'
distfile = prefix + suffix

data = np.load(distfile)
pocket_arr = data['pocket_arr']
pocket_XE1_dist_arr = data['pocket_XE1']
pocket_XE2_dist_arr = data['pocket_XE2']
pocket_XE3_dist_arr = data['pocket_XE3']
pocket_XE4_dist_arr = data['pocket_XE4']
pocket_XE5_dist_arr = data['pocket_XE5']
pocket_XE6_dist_arr = data['pocket_XE6']
pocket_XE7_dist_arr = data['pocket_XE7']
pocket_XE8_dist_arr = data['pocket_XE8']
pocket_XE9_dist_arr = data['pocket_XE9']

pockets = [pocket_XE1_dist_arr, pocket_XE2_dist_arr, pocket_XE3_dist_arr, pocket_XE4_dist_arr, pocket_XE5_dist_arr,
           pocket_XE6_dist_arr, pocket_XE7_dist_arr, pocket_XE8_dist_arr, pocket_XE9_dist_arr]


def smoothScan(pockets, frame, scan_dist, timesteps):
    active_pocket = []
    for i in range(scan_dist):
        if timesteps - 1 > frame + i:
            active_pocket.append(get_active_pocket(pockets, frame + i))
    if len(active_pocket) > 0:
        print(frame, mode(active_pocket), np.unique(active_pocket, return_counts=True))
        return mode(active_pocket)
    else:
        return 0


def get_active_pocket(pocketList, frame):
    compareList = []
    for pocket in pocketList:
        compareList.append(pocket[frame])
    if min(compareList) < 2.5:
        return compareList.index(min(compareList)) + 1
    else:
        return np.nan

num_timesteps = len(pocket_arr)

max_val = np.max([
    np.max(pocket_XE1_dist_arr), np.max(pocket_XE2_dist_arr), np.max(pocket_XE3_dist_arr), np.max(pocket_XE4_dist_arr),
    np.max(pocket_XE5_dist_arr), np.max(pocket_XE6_dist_arr), np.max(pocket_XE7_dist_arr)])

max_val = min(max_val, 20)

x = np.arange(num_timesteps)

# plot in which pocket the H2 is
nplot = 100

x = np.arange(num_timesteps)

SMALL_SIZE = 15
MEDIUM_SIZE = 22
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(1, 1, 1)
x = x / 1000

# # plot in which pocket the H2 is
# pocket_arr = []
# scan_dist = 2500
# for i in range(0, num_timesteps, scan_dist):
#     next100Pockets = smoothScan(pockets, i, scan_dist, num_timesteps)
#     pocket_arr += [next100Pockets]*scan_dist
#
# plt.plot(x+1000, pocket_arr, 'b', linewidth=2.0, label='pocket state')
# plt.plot(x[::nplot], pocket_arr[::nplot], 'ob', ms=3.0, linewidth=2.0)
#
# plt.ylabel("Pocket index")
# tbox = TextArea(
#     'B',
#     textprops=dict(
#         color='k', fontsize=40, ha='center', va='center')
# )

# plot all distances to H2 for debugging
plt.plot(x[::nplot], pocket_XE1_dist_arr[::nplot], linewidth=1.0, label='dist to XE1')
plt.plot(x[::nplot], pocket_XE2_dist_arr[::nplot], linewidth=1.0, label='dist to XE2')
plt.plot(x[::nplot], pocket_XE3_dist_arr[::nplot], linewidth=1.0, label='dist to XE3')
plt.plot(x[::nplot], pocket_XE4_dist_arr[::nplot], linewidth=1.0, label='dist to XE4')
plt.plot(x[::nplot], pocket_XE5_dist_arr[::nplot], linewidth=1.0, label='dist to XE5')
plt.plot(x[::nplot], pocket_XE6_dist_arr[::nplot], color='black', linewidth=2.0, label='dist to XE6')
plt.plot(x[::nplot], pocket_XE7_dist_arr[::nplot], color='gray', linewidth=1.0, label='dist to XE7')
plt.plot(x[::nplot], pocket_XE8_dist_arr[::nplot], color='gold', linewidth=1.0, label='dist to XE8')
plt.plot(x[::nplot], pocket_XE9_dist_arr[::nplot], color='yellow', linewidth=1.0, label='dist to XE9')
plt.xlim(0,600)
plt.legend(fancybox=True, loc='upper right', title='Distance H$_2$\nto Pocket', framealpha=1)
plt.ylabel("Distance ($\mathrm{\AA}$)")
tbox = TextArea(
    'D',
    textprops=dict(
        color='k', fontsize=40, ha='center', va='center')
    )


anchored_tbox = AnchoredOffsetbox(
    loc="upper right", child=tbox, pad=0., frameon=False,
    bbox_to_anchor=(0.05, 0.97),
    bbox_transform=ax.transAxes, borderpad=0.)
ax.add_artist(anchored_tbox)

plt.xlabel("Time (ps)")
plt.xlim(0, 925)
plt.ylim(0, 25)
# plt.yticks(np.arange(1, 10))

plt.savefig(dpi=400, fname="MorseMDCM_B", )
plt.show()
