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
    angle_errors = [5 for _ in angles]
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
    recoil_num_counts_error = []

    # load recoil counts from peak_info.json
    with open('data/peak_info.json', 'r') as f:
        data = json.load(f)
    for angle in angles:
        filename = f'recoil{angle}.Chn'
        num_counts = data[filename]['num_counts']
        recoil_num_counts.append(num_counts)
        num_counts_error = data[filename]['num_counts_err']
        recoil_num_counts_error.append(num_counts_error)

    print("recoil num counts", recoil_num_counts)

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
    recoil_flux_error = [recoil_num_counts_error[i] / live_times[i] for i in range(len(angles))]

    fluxes_dict = {}
    for i in range(len(angles)):
        fluxes_dict[f'recoil{angles[i]}'] = {
            'flux': recoil_flux[i],
            'flux_error': recoil_flux_error[i]
        }
    with open('data/fluxes.json', 'w') as f:
        json.dump(fluxes_dict, f, indent=4)

    # fit the flux vs angle to the Klein-Nishina formula
    popt, pcov = curve_fit(kn_function, angles, recoil_flux)
    print("KN popt", popt)
    # fit the flux vs angle to the Thomson formula
    popt_thomson, pcov_thomson = curve_fit(thomson, angles, recoil_flux)
    print("Thomson popt", popt_thomson)

    # add a 10% error to the fluxes
    recoil_flux_error = [0.1 * flux for flux in recoil_flux]

    # get chi squared values for the fits
    chi_squared_kn = np.sum((np.array(recoil_flux) - np.array(kn_function(angles, *popt)))**2 / np.array(recoil_flux_error)**2)
    chi_squared_thomson = np.sum((np.array(recoil_flux) - np.array(thomson(angles, *popt_thomson)))**2 / np.array(recoil_flux_error)**2)
    print(f"Klein-Nishina: chi^2 / dof = {chi_squared_kn} / {len(angles) - len(popt)}")
    print(f"Thomson: chi^2 / dof = {chi_squared_thomson} / {len(angles) - len(popt_thomson)}")

    # plot the fit
    plt.errorbar(angles, recoil_flux, xerr=angle_errors, yerr=recoil_flux_error, fmt='o', label='Fluxes with errors', color='blue')
    theta_plt = np.linspace(min(angles), max(angles), 100)
    plt.plot(theta_plt, kn_function(theta_plt, *popt), label=f"Klein-Nishina", color='blue')
    plt.plot(theta_plt, thomson(theta_plt, *popt_thomson), label=f"Thomson", color='red')
    plt.xlabel("Angle (degrees)", fontsize=14)
    plt.ylabel("Flux", fontsize=14)

    # annotate the chi squared values
    plt.text(0.4, 0.7, f"Klein-Nishina: $\chi^2$ / dof = {chi_squared_kn:.1f} / {len(angles) - len(popt)}", fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.4, 0.65, f"Thomson: $\chi^2$ / dof = {chi_squared_thomson:.1f} / {len(angles) - len(popt_thomson)}", fontsize=12, transform=plt.gca().transAxes)

    plt.title("Flux vs Angle", fontsize=16)
    plt.legend()
    plt.savefig("images/for_paper/klein_nishina.png")