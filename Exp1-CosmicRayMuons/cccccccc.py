'''
c stands for the speed of light
'''

import numpy as np
import matplotlib.pyplot as plt
from velocity_analysis import get_tof
import lin_regress as lr

distances = [100, 150, 200, 290]
corrected_distances = [112, 159, 207, 295]
distance_std = [9, 8, 6, 5]
# corrected_distances = [150, 200, 290]
file_names = [str(d) + 'cm_final' for d in distances]

# get the time of flight for each distance
tofs = [get_tof(file_name)[0] for file_name in file_names]
tof_errors = [get_tof(file_name)[1] for file_name in file_names]
# print(tofs)

# plot the time of flight vs distance
plt.plot(corrected_distances, tofs, 'o')
plt.xlabel('distance (cm)')
plt.ylabel('time of flight (ns)')
plt.title('time of flight vs distance')

# fit a line to the data

m, b, merr, berr = lr.lin_regression_with_error(corrected_distances, tofs,distance_std,  tof_errors, 1000)
max_line = [(m-merr)*d+(b+berr) for d in corrected_distances]
best_line = [(m)*d+(b) for d in corrected_distances]
min_line = [(m+merr)*d+(b-berr) for d in corrected_distances]



plt.errorbar(corrected_distances, tofs, xerr=distance_std, yerr=tof_errors,ls='none',marker="o")
plt.fill_between(corrected_distances, min_line, max_line,alpha=0.5)
plt.plot(corrected_distances, best_line,"b")


#plt.plot(distances, m*np.array(distances) + b)

plt.savefig('./images/tof_vs_distance.png')

# calculate the speed of light
# m is in units of ns/cm
c = (1/m) * 1e9 # in units of cm/s
cerr = (merr/m**2)*1e9
cc = 2.998e10 # in units of cm/s
print(f"fraction of speed of light: {c/cc} +- {cerr/cc}")

