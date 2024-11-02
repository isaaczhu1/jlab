'''

We will count the flux of the scattered photons at different angles,
and try to verify the Klein-Nishina formula for the differential cross section of Compton scattering.

'''

import numpy as np
import matplotlib.pyplot as plt

from load_data import *
from get_calib import *
from energy_angle import *
import json


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
    scatter_fluxes = [get_energy(f'goode/scatter{angle}.Chn', expected, verbose=False, return_flux=True) for angle, expected in zip(angles, expected_scatter)]
    recoil_fluxes = [get_energy(f'goode/recoil{angle}.Chn', expected, verbose=False, return_flux=True) for angle, expected in zip(angles, expected_recoil)]

    with open('data/live_times.json', 'r') as f:
        live_times = json.load(f)
        scatter_runtimes = [live_times[f'goode/scatter{angle}.Chn'] for angle in angles]
        recoil_runtimes = [live_times[f'goode/recoil{angle}.Chn'] for angle in angles]

    scatter_fluxes = [flux / runtime for flux, runtime in zip(scatter_fluxes, scatter_runtimes)]
    recoil_fluxes = [flux / runtime for flux, runtime in zip(recoil_fluxes, recoil_runtimes)]

    print("scatter fluxes:", scatter_fluxes)
    print("recoil fluxes:", recoil_fluxes)

    avg_fluxes = [r for s, r in zip(scatter_fluxes, recoil_fluxes)]
    print("average fluxes:", avg_fluxes)

    # plot the fluxes
    plt.scatter(angles, avg_fluxes, label="Scattered photons")

    plt.xlabel("Angle of detector")
    plt.ylabel("Flux (counts/s)")
    plt.legend()
    plt.savefig("images/for_paper/flux_angle.png")
