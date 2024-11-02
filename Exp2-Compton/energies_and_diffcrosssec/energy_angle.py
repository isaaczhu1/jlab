'''

We have compton scattering data at several angles.
The data is in the form of a histogram, with the x-axis being the energy of the scattered photon and the y-axis being the number of counts.
There is noise. We will only look at windows around the peaks.
Then we will fit the peaks to a gaussian distribution to find the peak energy.
Finally, we will plot the peak energy as a function of the angle of the detector.

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from load_data import *
from get_calib import *
import json
import pickle

def get_energy(filename, exp_eng_peak, verbose=False, name="none"):
    '''
    Given a file and an expected energy of peak, returns the energy of the peak
    '''
    if verbose:
        print("Looking at file", filename)
        print("Expected energy of peak:", exp_eng_peak)

    counts = read_data(filename)['counts']


    # convert the energy estimate to an index estimate
    with open('data/calib_x.pkl', 'rb') as f:
        data = pickle.load(f)
        counts_calib_x = data[filename]
    exp_int_peak = inv_calib_x(orig_x=range(len(counts)), calib_x=counts_calib_x, val=exp_eng_peak)
    if verbose:
        print(f"Expected energy of {exp_eng_peak} corresponds to index {exp_int_peak}")


    peak_info = get_peak(counts, exp_int_peak, debug=verbose, plot=True, peak_name=name)
    peak = peak_info['mean']
    peak_std = abs(peak_info['std'])
    peak_amplitude = peak_info['amplitude']

    # save the peak info to the json file data/peak_info.json
    with open('data/peak_info.json', 'r') as f:
        data = json.load(f)
    data[filename] = {
        "peak": peak,
        "std": peak_std,
        "amplitude": peak_amplitude,
        "num_counts": peak_std * peak_amplitude
    }
    with open('data/peak_info.json', 'w') as f:
        json.dump(data, f, indent=4)
        

    if verbose:
        print(f"The actual peak index: {peak}")

    energy = counts_calib_x[peak]

    return energy

if __name__ == "__main__":
    angles = [30, 60, 90, 120, 135]

    INCOMING_ENERGY = 661.66
    expected_scatter = [
        INCOMING_ENERGY / (1 + (INCOMING_ENERGY / 511) * (1 - np.cos(np.radians(angle))))
        for angle in angles
    ]
    expected_recoil = [
        INCOMING_ENERGY - energy
        for energy in expected_scatter
    ]
    
    verbose = True
    
    # scatter_energies = [get_energy(f'scatter{angle}.Chn', expected, verbose=verbose, name="scatter"+str(angle)) for angle, expected in zip(angles, expected_scatter)]
    # print("espected scattering energies:", expected_scatter)
    # print("actual scattering energies:", scatter_energies)

    # recoil_energies = [get_energy(f'recoil{angle}.Chn', expected, verbose=verbose, name="recoil"+str(angle)) for angle, expected in zip(angles, expected_recoil)]
    # print("expected recoil energies:", expected_recoil)
    # print("actual recoil energies:", recoil_energies)

    scatter_energies = []
    recoil_energies = []

    for angle in angles:
        filename = f'scatter{angle}.Chn'
        counts = read_data(filename)['counts']
        # smooth counts by replacing every bin with the average of the bin and its neighbors
        counts = np.array([np.mean(counts[max(0, i-1):min(len(counts), i+2)]) for i in range(len(counts))])
        scatter_energy = get_energy(f'scatter{angle}.Chn', expected_scatter[angles.index(angle)], verbose=verbose, name="scatter"+str(angle))
        print(f"Expected scattering energy: {expected_scatter[angles.index(angle)]}")
        print(f"Actual scattering energy: {scatter_energy}")
        plt.plot(range(len(counts)), counts, alpha=0.5)
        plt.title(f'Calibrated {filename}')
        plt.xlabel('Energy (keV)', fontsize=13)
        plt.ylabel('Counts', fontsize=13)
        plt.savefig(f'./images/energy_measurements/{filename[:-4]}_calib.png')
        plt.clf()
        scatter_energies.append(scatter_energy)

    for angle in angles:
        filename = f'recoil{angle}.Chn'
        counts = read_data(filename)['counts']
        # smooth counts by replacing every bin with the average of the bin and its neighbors
        counts = np.array([np.mean(counts[max(0, i-1):min(len(counts), i+2)]) for i in range(len(counts))])
        recoil_energy = get_energy(f'recoil{angle}.Chn', expected_recoil[angles.index(angle)], verbose=verbose, name="recoil"+str(angle))
        print(f"Expected recoil energy: {expected_recoil[angles.index(angle)]}")
        print(f"Actual recoil energy: {recoil_energy}")
        plt.plot(range(len(counts)), counts, alpha=0.5)
        plt.title(f'Calibrated {filename}')
        plt.xlabel('Energy (keV)', fontsize=13)
        plt.ylabel('Counts', fontsize=13)
        plt.savefig(f'./images/energy_measurements/{filename[:-4]}_calib.png')
        plt.clf()
        recoil_energies.append(recoil_energy)

    data = {
        "angles": angles,
        "expected_scatter": expected_scatter,
        "expected_recoil": expected_recoil,
        "scatter_energies": scatter_energies,
        "recoil_energies": recoil_energies
    }

    with open('data/energy_angle.json', 'w') as f:
        json.dump(data, f, indent=4)


    