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

def get_energy(filename, exp_eng_peak, verbose=False, return_flux=False, name="none"):
    '''
    Given a file and an expected energy of peak, returns the energy of the peak
    '''
    if verbose:
        print("Looking at file", filename)
        print("Expected energy of peak:", exp_eng_peak)

    counts = read_data(filename)['counts']


    # convert the energy estimate to an index estimate
    counts_calib_x = get_calib_x(filename)
    exp_int_peak = inv_calib_x(orig_x=range(len(counts)), calib_x=counts_calib_x, val=exp_eng_peak)
    if verbose:
        print(f"Expected energy of {exp_eng_peak} corresponds to index {exp_int_peak}")


    peak = get_peak(counts, exp_int_peak, debug=verbose, return_flux=return_flux, name=name)

    if return_flux:
        return peak

    if verbose:
        print(f"The actual peak index: {peak}")

    energy = counts_calib_x[peak]

    return energy

if __name__ == "__main__":
    angles = [30, 60, 90, 120]
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
    scatter_energies = [get_energy(f'goode/scatter{angle}.Chn', expected, verbose=verbose, name="scatter"+str(angle)) for angle, expected in zip(angles, expected_scatter)]
    print("espected scattering energies:", expected_scatter)
    print("actual scattering energies:", scatter_energies)

    recoil_energies = [get_energy(f'goode/recoil{angle}.Chn', expected, verbose=verbose, name="recoil"+str(angle)) for angle, expected in zip(angles, expected_recoil)]
    print("expected recoil energies:", expected_recoil)
    print("actual recoil energies:", recoil_energies)

    data = {
        "angles": angles,
        "expected_scatter": expected_scatter,
        "scatter_energies": scatter_energies,
        "expected_recoil": expected_recoil,
        "recoil_energies": recoil_energies
    }

    with open('data/energy_angle.json', 'w') as f:
        json.dump(data, f, indent=4)


    