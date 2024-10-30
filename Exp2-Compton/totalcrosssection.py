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
import math
import random

rebin_amount = 2 # How much rebinning is done.

## Load data files

total_counts = [0,0]
total_errors = [0,0]
live_times = []
attenuation = [0,2.48,4.98,7.5,10.03,12.50]
for j in range(2, 4):
    data_filenames = [f"cs{j}_{i}" for i in range(6)]
    #background_filename = "cs3_background"
    
    
    data = [read_data(filename) for filename in data_filenames] 
    #data.append(read_data(background_filename))
    live_times = [data[i]['time'] for i in range(6)]
    
    #print(live_times)
    
    raw_counts = np.array([rebin(data[i]['counts'],rebin_amount) for i in range(6)])
    #background = np.array(rebin(data[6]['counts'],rebin_amount))
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
    
    
    print(window_start)
    print(window_end)
    total_counts[j-2] = [sum(raw_counts[i][window_start:window_end]) for i in range(6)]
    total_errors[j-2] = np.sqrt([sum(raw_count_errors[i][window_start:window_end]**2) for i in range(6)])
#plt.plot([window_start, window_end],[0,0],label="window")

#plt.legend()
## Collect the count rates

"""
total_counts = [sum(raw_counts[i]) for i in range(6)]
total_errors = np.sqrt([sum(raw_count_errors[i]**2) for i in range(6)])
"""
final_counts = []
final_errors = []
final_attenuation = []
for i in range(6):
    final_counts.append(total_counts[0][i]/total_counts[0][0])
    final_counts.append(total_counts[1][i]/total_counts[1][0])
    final_errors.append(total_errors[0][i]/total_counts[0][0])
    final_errors.append(total_errors[1][i]/total_counts[1][0])
    final_attenuation.append(attenuation[i])
    final_attenuation.append(attenuation[i])

plt.errorbar(final_attenuation, final_counts,xerr=0.01,yerr=final_errors,ls="None",marker=".")
## Fit an exponential

exp = lambda x, A, B, C: A*np.exp(-B*x)+C

guess = [total_counts[0], 0.08, 150]
        
[A, mu, C], cov = scipy.optimize.curve_fit(exp, final_attenuation, final_counts, p0=guess, sigma=final_errors)

plt.plot(attenuation, [exp(x, A, mu, C) for x in attenuation])

## Data analysis junk

print(f"Decay coefficent: {mu}+- {np.sqrt(cov[1][1])}")
print(f"Background: {C}+- {np.sqrt(cov[2][2])}")

print(f"Cross-section: {mu/(3.36*10**(27))} m^2")

chi2 = sum([((total_counts[i] - exp(attenuation[i], A, mu,C))/(total_errors[i]))**2 for i in range(6)])
print(f"Chi-squared: {chi2} / 4")

