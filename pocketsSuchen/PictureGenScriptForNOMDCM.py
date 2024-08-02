import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

distfile = 'pocket_arrays.npz'


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


max_val = np.max([
    np.max(pocket_XE1_dist_arr), np.max(pocket_XE2_dist_arr), np.max(pocket_XE3_dist_arr), np.max(pocket_XE4_dist_arr),
    np.max(pocket_XE5_dist_arr), np.max(pocket_XE6_dist_arr), np.max(pocket_XE7_dist_arr)])

num_timesteps = 10000


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


# plot in which pocket the H2 is
# plt.plot(x, pocket_arr, linewidth=2.0, label='pocket state')
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(1, 1, 1)
x = x/10





# plt.plot(x, pocket_arr, linewidth=2.0)
# plt.ylabel("Pocket index")
# tbox = TextArea(
#     'A',
#     textprops=dict(
#         color='k', fontsize=40, ha='center', va='center')
#     )


# plot all distances to H2 for debugging
plt.plot(x, pocket_XE1_dist_arr, linewidth=1.0, label=' 1')
plt.plot(x, pocket_XE2_dist_arr, linewidth=1.0, label=' 2')
plt.plot(x, pocket_XE3_dist_arr, linewidth=1.0, label=' 3')
plt.plot(x, pocket_XE4_dist_arr, linewidth=1.0, label=' 4')
plt.plot(x, pocket_XE5_dist_arr, linewidth=1.0, label=' 5')
plt.plot(x, pocket_XE6_dist_arr, color='black', linewidth=2.0, label=' 6')
plt.plot(x, pocket_XE7_dist_arr, color='gray', linewidth=1.0, label=' 7')
plt.plot(x, pocket_XE8_dist_arr, color='gold', linewidth=1.0, label=' 8')
plt.plot(x, pocket_XE9_dist_arr, color='yellow', linewidth=1.0, label=' 9')
plt.ylabel("Distance ($\mathrm{\AA}$)")
plt.legend(fancybox=True, loc='upper left', title='Distance H$_2$\nto Pocket', framealpha=1)
tbox = TextArea(
    'C',
    textprops=dict(
        color='k', fontsize=40, ha='center', va='center')
    )



anchored_tbox = AnchoredOffsetbox(
    loc="upper right", child=tbox, pad=0., frameon=False,
    bbox_to_anchor=(0.97, 0.97),
    bbox_transform=ax.transAxes, borderpad=0.)
ax.add_artist(anchored_tbox)


plt.xlim(0, 600)
plt.xlabel("Time (ps)")

plt.savefig(dpi=400, fname="CGenFF", )
plt.show()

