import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from load_data import *
from get_calib import *
import json
import pickle

def generate_background():
    angles = [30, 60, 90, 120, 135]

    all_counts = {}

    for angle in angles:
        filename = f"recoil{angle}.Chn"
        counts = read_data(filename)['counts']
        # smooth counts by replacing every bin with the average of the bin and its neighbors
        N = 2
        counts = np.array([np.mean(counts[max(0, i+1-N):min(len(counts), i+N)]) for i in range(len(counts))])
        with open('data/calib_x.pkl', 'rb') as f:
            data = pickle.load(f)
            counts_calib_x = data[f'recoil{angle}.Chn']


        # print(len(counts))

        # find the average counts between 0.62*len(counts) 0.78*len(counts)
        L = int(0.62*len(counts))
        R = int(0.78*len(counts))
        avg = np.mean(counts[L:R])
        # normalize the data so avg = 1
        counts = counts / avg

        all_counts[angle] = counts

        plt.plot(counts_calib_x, counts, alpha=0.5, label=f'Recoil at {angle} degrees')
        # print(max(counts_calib_x))

    # at every index, find the average of the two smallest counts
    background = []
    for i in range(len(counts)):
        vals = [all_counts[angle][i] for angle in angles]
        vals.sort()
        if i < 0.62*len(counts):
            vals = vals[:3]
        background.append(np.mean(vals))

    background = np.array(background)
    plt.plot(counts_calib_x, background, alpha=1, color='black', label='Background')

    plt.xlabel('Energy (keV)', fontsize=14)
    plt.ylabel('Normalized counts', fontsize=14)
    plt.title('Measuring the Recoil background', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()

    plt.savefig('images/recoil_alignment/getting_background.png')
    plt.show()

    # pickle the background
    # with open('data/background.pkl', 'wb') as f:
    #     pickle.dump(background, f)

generate_background()


with open('data/background.pkl', 'rb') as f:
    background = pickle.load(f)

angles = [30, 60, 90, 120, 135]
new_counts = {}

# now, plot the background subtracted data
for angle in angles:
    filename = f"recoil{angle}.Chn"
    counts = read_data(filename)['counts']
    # smooth counts by replacing every bin with the average of the bin and its neighbors
    N = 2
    counts = np.array([np.mean(counts[max(0, i+1-N):min(len(counts), i+N)]) for i in range(len(counts))])
    with open('data/calib_x.pkl', 'rb') as f:
        data = pickle.load(f)
        counts_calib_x = data[f'recoil{angle}.Chn']

    # find the average counts between 0.62*len(counts) 0.78*len(counts)
    L = int(0.62*len(counts))
    R = int(0.78*len(counts))
    avg = np.mean(counts[L:R])
    # normalize the data so avg = 1
    counts = counts / avg

    # subtract the background
    counts = counts - background

    new_counts[angle] = counts

    PLT_ANGLE = 90
    if True: #angle == PLT_ANGLE:
        # fit a gaussian to (counts_calib_x, counts) around the peak
        peak_index = np.argmax(counts)
        peak_val = counts[peak_index]
        def gaussian(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

        popt_guess = [peak_val, counts_calib_x[peak_index], len(counts) / 10]
        popt, _ = curve_fit(gaussian, counts_calib_x, counts, p0=popt_guess)

        a, x0, sigma = popt
        sigma = abs(sigma)
        # scale up the amplitude to match the original counts
        a = a * avg
        a_err = np.sqrt(np.diag(_)[0])
        x0_err = np.sqrt(np.diag(_)[1])
        sigma_err = np.sqrt(np.diag(_)[2])
        peak_info = {
            "mean": x0,
            "std": sigma,
            "amplitude": a,
            "num_counts": a * sigma, # up to a constant factor
            "num_counts_err": np.sqrt((a_err * sigma)**2 + (a * sigma_err)**2)
        }

        # dump into new_peak_info.json
        with open('data/new_peak_info.json', 'r') as f:
            data = json.load(f)
        data[filename] = peak_info
        with open('data/new_peak_info.json', 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Stored peak info for {filename} in data/new_peak_info.json")

        
        plt.plot(counts_calib_x, counts, alpha=0.5, label=f'Recoil at {angle} degrees')
        plt.plot(counts_calib_x, gaussian(counts_calib_x, *popt), label=f'Gaussian fit around peak', color='red', linewidth=2)
    # print(max(counts_calib_x))

plt.xlabel('Energy (keV)', fontsize=14)
plt.ylabel('Normalized counts', fontsize=14)
# plt.title('Background subtracted data', fontsize=16)
plt.title(f'Background subtracted data: Recoil at {PLT_ANGLE} degrees', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.legend()
# plt.savefig(f'images/recoil_alignment/fitting_subtracted_{PLT_ANGLE}.png')
plt.show()

assert False

# find the peaks in the new counts


peaks = {}


INCOMING_ENERGY = 661.66

expected_scatter = [
        INCOMING_ENERGY / (1 + (INCOMING_ENERGY / 511) * (1 - np.cos(np.radians(angle))))
        for angle in angles
    ]

expected_recoil = [
        INCOMING_ENERGY - energy
        for energy in expected_scatter
    ]

for i, angle in enumerate(angles):
    # print("hi", angle)
    expected_energy = expected_recoil[i]
    counts = new_counts[angle]
    with open('data/calib_x.pkl', 'rb') as f:
        data = pickle.load(f)
        counts_calib_x = data[filename]
    exp_int_peak = inv_calib_x(orig_x=range(len(counts)), calib_x=counts_calib_x, val=expected_energy)


    # print("Expected index of peak:", exp_int_peak)
    print("Expected energy of peak:", expected_energy)

    # plt.show()
    plt.plot(counts)

    # find the peak
    peak_info = get_peak(counts, exp_int_peak, debug=False, plot=True, peak_name=f"bonk{angle}")
    # print("actual peak:" , peak_info["mean"])
    print("actual peak energy:", counts_calib_x[peak_info["mean"]])
# plt.show()