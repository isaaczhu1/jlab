'''

load the data in data/energy_angle.json and plot the expected vs actual energies for scatter and recoil

'''

import json
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    with open("data/energy_angle.json", "r") as f:
        data = json.load(f)

    angles = data['angles']
    expected_scatter = data['expected_scatter']
    scatter_energies = data['scatter_energies']
    expected_recoil = data['expected_recoil']
    recoil_energies = data['recoil_energies']

    plt.scatter(angles, expected_scatter, label="Expected scatter energy", color='blue')
    plt.scatter(angles, scatter_energies, label="Actual scatter energy", color='blue', marker='x')
    plt.scatter(angles, expected_recoil, label="Expected recoil energy", color='green')
    plt.scatter(angles, recoil_energies, label="Actual recoil energy", color='green', marker='x')

    # plot the sum of the energies
    plt.scatter(angles, [expected_scatter[i] + expected_recoil[i] for i in range(len(angles))], label="Expected sum of energies", color='red')
    plt.scatter(angles, [scatter_energies[i] + recoil_energies[i] for i in range(len(angles))], label="Actual sum of energies", color='red', marker='x')

    plt.xlabel("Angle of detector")
    plt.ylabel("Energy (keV)")
    plt.legend()
    plt.savefig("images/for_paper/energy_angle.png")

    plt.clf()

    # plot the inverse energy vs 1 - cos(theta)
    plt.scatter([1 - np.cos(np.radians(angle)) for angle in angles], [1 / energy for energy in scatter_energies], label="Scattered photons")

    # draw a best fit line
    m, b = np.polyfit([1 - np.cos(np.radians(angle)) for angle in angles], [1 / energy for energy in scatter_energies], 1)

    x = np.linspace(0, max([1 - np.cos(np.radians(angle)) for angle in angles]), 100)
    y = m * x + b

    plt.plot(x, y, label=f"y = {m:.4f}x + {b:.4f}")

    plt.xlabel("1 - $\\cos(\\theta)$ (dimensionless)")
    plt.ylabel("Inverse energy (1/keV)")
    plt.ylim(0, max([1 / energy for energy in scatter_energies])*1.1)
    plt.legend()
    plt.savefig("images/for_paper/inverse_energy_angle.png")

    plt.clf()