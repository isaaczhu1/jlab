'''
c stands for the speed of light
'''

import numpy as np
import matplotlib.pyplot as plt
from velocity_analysis import get_tof
import lin_regress as lr
from distance_distribution_MC import get_mean_dist
import pickle

distances = [100, 150, 200, 250, 290]



# if the pickled file exists, load it, otherwise run the simulation
try:
    with open('./corrected_distances.pkl', 'rb') as f:
        corrected_distances = pickle.load(f)
        corrected_distances_std = pickle.load(f)
    print(f"Loaded corrected distances from pickle: {corrected_distances}, {corrected_distances_std}")
except:
    # perform the Monte Carlo simulation to correct for the mean flight distance
    MC_SIM_TRIALS = 10000
    print("STARTING MC SIMULATION")
    sim = [get_mean_dist(d, MC_SIM_TRIALS, verbose=True) for d in distances]
    print("MC SIMULATION FINISHED")
    corrected_distances = [s[0] for s in sim]
    corrected_distances_std = [s[1] for s in sim]
    # pickle the corrected distances
    with open('./corrected_distances.pkl', 'wb') as f:
        pickle.dump(corrected_distances, f)
        pickle.dump(corrected_distances_std, f)


file_names = [str(d) + 'cm_final' for d in distances]

# get the time of flight for each distance
tofs = [get_tof(file_name)[0] for file_name in file_names]
tof_errors = [get_tof(file_name)[1] for file_name in file_names]

print(f"Loaded time-of-flights: {tofs}, {tof_errors}")
# print(tofs)

# plot the time of flight vs distance
# plt.plot(corrected_distances, tofs, 'o')
plt.xlabel('Distance (cm)')
plt.ylabel('Time of Flight (ns)')
plt.title('Time of Flight vs Distance')

# fit a line to the data

m, b, merr, berr = lr.lin_regression_with_error(corrected_distances, tofs, corrected_distances_std,  tof_errors, 1000)
max_line = [(m-merr)*d+(b+berr) for d in corrected_distances]
best_line = [(m)*d+(b) for d in corrected_distances]
min_line = [(m+merr)*d+(b-berr) for d in corrected_distances]



plt.errorbar(corrected_distances, tofs, xerr=corrected_distances_std, yerr=tof_errors, ls='none', marker=".", color = 'tab:orange',label="Data")
plt.fill_between(corrected_distances, min_line, max_line,alpha=0.3,label="$1\sigma$ error")
plt.plot(corrected_distances, best_line, color = 'tab:blue',label="Best fit")
plt.legend(loc="lower right")


# calculate chi-squared of trend line
chi = sum([(tofs[i]-best_line[i])**2/(tof_errors[i]**2) for i in range(len(tofs))])
print(chi)
# calculate the speed of light
# m is in units of ns/cm
c = (1/m) * 1e7 # in units of m/s
cerr = (merr/m**2)*1e7
cc = 2.998e8 # in units of m/s
print(f"fraction of speed of light: {c/cc} +- {cerr/cc}")

# divide c by 1e8 to make it easier to read
c_text = f"({c/1e8:.2f} $\pm$ {cerr/1e8:.2f}) x $10^8$ m/s"
plt.text(110, 35, f"Measured velocity of muons = {c_text}")
plt.text(110, 34, f"Fraction of c: {c/cc:.2f} $\pm$ {cerr/cc:.2f}")
plt.text(110, 33, f"$\chi^2$: {chi:.2f} / 3, probability: 12%")


plt.savefig('./images/tof_vs_distance.png')

