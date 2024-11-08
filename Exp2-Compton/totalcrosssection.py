# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:46:25 2024

@author: evanb
"""

from load_data import read_data, rebin
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
import math
import random

rebin_amount = 2 # How much rebinning is done.

## Load data files

total_counts = [0,0]
total_errors = [0,0]
live_times = []
attenuation = [0,2.48,4.98,7.5,10.03,12.50]
data_filenames = [f"attenuation/cs2_{i}.Chn" for i in range(6)]
background_filename = "attenuation/cs2_background.Chn"
    
data = [read_data(filename) for filename in data_filenames] 
data.append(read_data(background_filename))
live_times = [data[i]['time'] for i in range(6)]

#print(live_times)

raw_counts = np.array([rebin(data[i]['counts'],rebin_amount) for i in range(6)])
background = np.array(rebin(data[6]['counts'],rebin_amount))
"""
deadspace = 0 # Index where the CFD cutoff elminates data
while raw_counts[0][deadspace] == 0:
    deadspace += 1

for i in range(6):
    for j in range(deadspace):
        raw_counts[i][j] = raw_counts[i][j+deadspace]
        
## Methodology above is to fill the deadspace to the CFD with that level of signal.
"""
raw_count_errors = np.sqrt(raw_counts)
#background_error = np.sqrt(background)

#background /= live_times[6]
#background_error /= live_times[6]

#plt.plot(raw_counts[3])

for i in range(6):
    #plt.plot(raw_counts[i],label=i)
    raw_counts[i] /= live_times[i]
    raw_count_errors[i] /= live_times[i]
    


## Create a window to analyize the data. Methodology is when there is no longer 
# a singificant distinction (3 sigma) between the no and max attenuation

window_start = 1550
window_end = 1550

sigma = (raw_counts[0][window_start]-raw_counts[-1][window_start])/(np.sqrt(raw_count_errors[0][window_start]**2 + raw_count_errors[-1][window_start]**2))
while sigma > 3:
    window_start -= 1
    sigma = (raw_counts[0][window_start]-raw_counts[-1][window_start])/(np.sqrt(raw_count_errors[0][window_start]**2 + raw_count_errors[-1][window_start]**2))

sigma = (raw_counts[0][window_end]-raw_counts[-1][window_end])/(np.sqrt(raw_count_errors[0][window_end]**2 + raw_count_errors[-1][window_end]**2))
while sigma > 3:
    window_end += 1
    sigma = (raw_counts[0][window_end]-raw_counts[-1][window_end])/(np.sqrt(raw_count_errors[0][window_end]**2 + raw_count_errors[-1][window_end]**2))


total_counts = [sum(raw_counts[i][window_start:window_end]) for i in range(6)]
total_errors = np.sqrt([sum(raw_count_errors[i][window_start:window_end]**2) for i in range(6)])


#total_counts = [sum(raw_counts[i]) for i in range(6)]
#total_errors = np.sqrt([sum(raw_count_errors[i]**2) for i in range(6)])

difference = []
for i in range(len(raw_counts[0])):
    if (np.sqrt(raw_count_errors[0][i]**2 + raw_count_errors[-1][i]**2)) > 0: 
        difference.append((raw_counts[0][i]-raw_counts[-1][i])/(np.sqrt(raw_count_errors[0][i]**2 + raw_count_errors[-1][i]**2)))
    else:
        difference.append(0)
fig, ax =plt.subplots()
ax.plot(difference)
#ax.plot(raw_counts[0],label="No Attenuation")
#ax.plot(raw_counts[5], label="Max Attenuation")
rect1 = plt.Rectangle((window_start, 0), (window_end-window_start), 25, fill=False, color='r',linewidth=2)
ax.add_patch(rect1)
plt.xlabel("Voltage Bin",size=14)
plt.ylabel("Deviation between max and min attenuation ($\sigma$)",size=14)
plt.title("Deviation and Selection Window")
#plt.plot([window_start, window_end],[0,0],label="window")

plt.legend()
#plt.savefig('./images/selection.png')
## Collect the count rates

"""
total_counts = [sum(raw_counts[i]) for i in range(6)]
total_errors = np.sqrt([sum(raw_count_errors[i]**2) for i in range(6)])



## Fit an exponential

exp = lambda x, A, B, C: A*np.exp(-B*x)+C

guess = [total_counts[0], 0.08, 0]
        
[A, mu, C], cov = scipy.optimize.curve_fit(exp, attenuation, total_counts, p0=guess, sigma=total_errors)

plt.errorbar(attenuation, [total_counts[i]-exp(attenuation[i], A, mu, C) for i in range(len(attenuation))],xerr=0.01,yerr=total_errors,ls="None",marker=".",label="Data Residuals")
plt.plot([-0.5,13],[0,0],ls="--")
chi2 = sum([((total_counts[i] - exp(attenuation[i], A, mu,C))/(total_errors[i]))**2 for i in range(6)])

#plt.plot(range(14), [exp(x, A, mu, C) for x in range(14)], label="Exponential Fit")
plt.legend()
plt.title("Residuals of Data")
plt.xlabel("Plastic Layer Thickness (cm)",size=12)
plt.ylabel("Residual (Hz)",size=12)
#plt.savefig("./images/Residual.png")

#plt.annotate(f"Decay Coefficient: {mu:.3f} $\pm$ {np.sqrt(cov[1][1]):.3f} 1/cm", (5.7,500))
#plt.annotate(f"Background Level: {C:.1f} $\pm$ {np.sqrt(cov[2][2]):.1f} Hz", (6.9,450))
#plt.annotate(f"$\chi^2$: {chi2:.2f}  /  3",(10.9,400))

#plt.savefig("./images/attenuation_curve.png",format="png")
## Data analysis junk
#print(total_counts, total_errors)
print(f"Cross-section: ({mu/(3.36*10**(-2)):.2f}+-{np.sqrt(cov[1][1])/(3.36*10**(-2)):.2f})e-29 m^2")

"""