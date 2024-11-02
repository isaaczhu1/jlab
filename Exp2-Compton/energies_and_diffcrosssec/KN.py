'''

We will count the flux of the scattered photons at different angles,
and try to verify the Klein-Nishina formula for the differential cross section of Compton scattering.

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from load_data import *
from get_calib import *
from energy_angle import *
import json

EPS = 661.66/511

def thomson(theta, r_e, c):
    return r_e**2 * (1 + np.cos(np.radians(theta))**2)/2 + c

def lam_ratio(theta):
    return 1/(1 + EPS*(1 - np.cos(np.radians(theta))))

def kn_function(theta, r_e, c):
    '''
    r_e: classical electron radius
    c: constant background
    '''
    lam = lam_ratio(theta)
    return 0.5 * r_e**2 * lam**2 \
        * (lam + 1/lam - np.sin(np.radians(theta))**2) + c


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

    recoil_num_counts = []

    # load recoil counts from peak_info.json
    with open('data/peak_info.json', 'r') as f:
        data = json.load(f)
    for angle in angles:
        filename = f'recoil{angle}.Chn'
        num_counts = data[filename]['num_counts']
        recoil_num_counts.append(num_counts)

    print(recoil_num_counts)

    # load live times from data/live_times.json
    live_times = []
    with open('data/live_times.json', 'r') as f:
        data = json.load(f)
    for angle in angles:
        filename = f'recoil{angle}.Chn'
        live_time = data[filename]
        live_times.append(live_time)

    # divide the number of counts by the live time to get the flux
    recoil_flux = [recoil_num_counts[i] / live_times[i] for i in range(len(angles))]

    fluxes_dict = dict(zip(angles, recoil_flux))
    with open('data/fluxes.json', 'w') as f:
        json.dump(fluxes_dict, f, indent=4)

    # fit the flux vs angle to the Klein-Nishina formula
    popt, pcov = curve_fit(kn_function, angles, recoil_flux)
    print(popt)
    # fit the flux vs angle to the Thomson formula
    popt_thomson, pcov_thomson = curve_fit(thomson, angles, recoil_flux)
    print(popt_thomson)
    # plot the fit
    plt.scatter(angles, recoil_flux)
    theta_plt = np.linspace(min(angles), max(angles), 100)
    plt.plot(theta_plt, kn_function(theta_plt, *popt), label=f"Klein-Nishina")
    plt.plot(theta_plt, thomson(theta_plt, *popt_thomson), label=f"Thomson")
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Flux")
    plt.title("Flux vs Angle")
    plt.legend()
    plt.savefig("images/for_paper/klein_nishina.png")