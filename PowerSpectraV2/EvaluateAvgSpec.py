import numpy as np

# Morse MDCM
# x_lim_min = 4407.
# x_lim_max = 4527.

# CGenFF
# x_lim_min = 4013.
# x_lim_max = 4133.


load_from = "no"
if load_from == "MDCMmorse":
    loader = np.load("spec_avg_saveMDCMmorse.npz")

    avg_spec = loader["spec_avg"]
    freq = loader["freq"]
    pocket = []
    freq2 = freq[np.logical_and(freq > 4407, freq < 4527)]

    middlepoint_list = []

    for i in range(len(avg_spec)):
        pocket.append(avg_spec[i][np.logical_and(freq > 4407, freq < 4527)])
        # print(pocket)
        middlepoint_list.append([freq2[pocket[i].argmax()], np.sum(pocket[i]*freq2)/np.sum(pocket[i])])
        print(freq2[pocket[i].argmax()])
        print(np.sum(pocket[i]*freq2)/np.sum(pocket[i]))

    np.savez("middle_mdcm_morse", pos=np.array(middlepoint_list))

else:
    loader = np.load("spec_avg_saveNOMDCM.npz")

    avg_spec = loader["spec_avg"]
    freq = loader["freq"]
    pocket = []
    freq2 = freq[np.logical_and(freq > 4013, freq < 4133)]

    middlepoint_list = []

    for i in range(len(avg_spec)):
        pocket.append(avg_spec[i][np.logical_and(freq > 4013, freq < 4133)])
        # print(pocket)
        middlepoint_list.append([freq2[pocket[i].argmax()], np.sum(pocket[i] * freq2) / np.sum(pocket[i])])
        print(freq2[pocket[i].argmax()])
        print(np.sum(pocket[i] * freq2) / np.sum(pocket[i]))

    np.savez("middle_NoMDCM", pos=np.array(middlepoint_list))