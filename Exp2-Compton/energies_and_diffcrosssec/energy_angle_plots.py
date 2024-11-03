'''

load the data in data/energy_angle.json and plot the expected vs actual energies for scatter and recoil

'''

import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR

if __name__ == "__main__":
    with open("data/energy_angle.json", "r") as f:
        data = json.load(f)

    # FIRST PART: PLOT ENERGY VS ANGLE

    angles = data['angles']
    angle_errors = [5 for _ in angles]
    expected_scatter = data['expected_scatter']
    scatter_energies = data['scatter_energies']
    scatter_energies_errors = data['scatter_energies_errors']
    expected_recoil = data['expected_recoil']
    recoil_energies = data['recoil_energies']
    recoil_energies_errors = data['recoil_energies_errors']

    theta_plt = np.linspace(min(angles), max(angles), 100)
    scatter_expected_plt = [661.66 / (1 + (661.66 / 511) * (1 - np.cos(np.radians(angle)))) for angle in theta_plt]
    recoil_expected_plt = [661.66 - energy for energy in scatter_expected_plt]

    plt.plot(theta_plt, scatter_expected_plt, label="Expected scatter", color='blue', linestyle='--')
    plt.errorbar(angles, scatter_energies, xerr=angle_errors, yerr=scatter_energies_errors, fmt='s', label="Actual scatter", color='blue')
    plt.plot(theta_plt, recoil_expected_plt, label="Expected recoil", color='green', linestyle='--')
    plt.errorbar(angles, recoil_energies, xerr=angle_errors, yerr=recoil_energies_errors, fmt='s', label="Actual recoil", color='green')

    # plot the sum of the energies
    plt.axhline(y=661.66, color='red', linestyle='--', label="Expected sum (661.66 keV)")
    sum_energies = [scatter + recoil for scatter, recoil in zip(scatter_energies, recoil_energies)]
    sum_energies_errors = [np.sqrt(scatter_err**2 + recoil_err**2) for scatter_err, recoil_err in zip(scatter_energies_errors, recoil_energies_errors)]
    plt.errorbar(angles, sum_energies, yerr=sum_energies_errors, fmt='s', label="Actual sum", color='red')

    sums = np.array(sum_energies)
    # get mean and standard deviation of the sum of the energies
    mean = np.mean(sums)
    std = np.std(sums)
    print(f"Mean of sum of energies: {mean:.1f} ± {std:.1f} keV, fractional error: {std/mean:.3f}")

    plt.xlabel("Angle of detector (degrees)", fontsize=14)
    plt.ylabel("Energy (keV)", fontsize=14)
    plt.title("Scattering/recoil energies vs angle", fontsize=16)

    # format the legend to lie above the line y=661.66 with two different columns
    plt.legend(loc='lower right', ncol=2)
    plt.savefig("images/for_paper/energy_angle.png")
    plt.clf()



    # SECOND PART: PLOT INVERSE SCATTERING ENERGY VS 1 - COS(THETA)

    # plot the inverse energy vs 1 - cos(theta)
    x_scatter = [1 - np.cos(np.radians(angle)) for angle in angles]
    x_scatter_errors = [np.sin(np.radians(angle)) * np.radians(5) for angle in angles]
    y_scatter = [1 / energy for energy in scatter_energies]
    y_scatter_errors = [scatter_energies_errors[i] / scatter_energies[i]**2 for i in range(len(scatter_energies))]

    plt.errorbar(
        [1 - np.cos(np.radians(angle)) for angle in angles],
        [1 / energy for energy in scatter_energies],
        xerr=x_scatter_errors,
        yerr=y_scatter_errors,
        fmt='o'
    )

    # draw a best fit line, using scipy odr since the errors are in both x and y
    def f(B, x):
        return B[0] * x + B[1]
    
    linear = Model(f)

    data = RealData(x_scatter, y_scatter, sx=x_scatter_errors, sy=y_scatter_errors)
    odr = ODR(data, linear, beta0=[1, 0])
    out = odr.run()
    m, c = out.beta
    m_err, c_err = out.sd_beta

    y = [1/energy for energy in scatter_energies]
    x_pred = [(y_val - c) / m for y_val in y]

    # get a chi squared value from the inverse plot
    chi_squared = sum([(x - x_pred[i])**2 / x_scatter_errors[i]**2 for i, x in enumerate(x_scatter)])
    # print(f"chi^2 / dof = {chi_squared:.1f} / {len(angles) - 2}")
    plt.text(0.1, 0.8, f"$\chi^2$ / dof = {chi_squared:.1f} / {len(angles) - 2}", transform=plt.gca().transAxes, fontsize=12)

    # the slope of the line is 1/m_e
    m_e = 1 / m
    m_e_err = m_e * m_err / m
    print(f"Obtained electron rest mass: {m_e:.1f} ± {m_e_err:.1f} keV")

    plt.plot(x_pred, y, label=f"Best fit line: $y = {m:.2f}x + {c:.2f}$")

    plt.title("Inverse scattering energy vs 1 - $\\cos(\\theta)$", fontsize=16)
    plt.xlabel("1 - $\\cos(\\theta)$ (dimensionless)", fontsize=14)
    plt.ylabel("Inverse energy (1/keV)", fontsize=14)
    plt.ylim(0, max([1 / energy for energy in scatter_energies])*1.1)
    plt.savefig("images/for_paper/inverse_energy_angle.png")
    plt.clf()