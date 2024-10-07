"""
NOTE: Most of this functionality can now be found in load_data.py

Created on Sun Sep 29 16:34:05 2024

@author: evanb, isaaczhu
"""


import numpy as np
import matplotlib.pyplot as plt

def read_data(filename):
    with open(f'./data/{filename}.Chn', 'rb') as f:
        data = f.read()
        # the data type is INTEGERS
        counts = np.frombuffer(data, dtype=np.uint16)

    print("counts.shape:",counts.shape)

    # must remove garbage stuff at the beginning and end
    counts = counts[16:-256]

    return counts


filename = 'calibration'
#preamble = read_data_preamble(filename)
counts = read_data(filename)

t0 = 0 #Starting bin time
tstep = 1 # 50/len(counts)  # Bin size
cutoff = 20 # Number of counts required to be considered in max

spacing = 20 # int(np.floor(1/tstep)) # Must be under the spacing between peaks in bins.

counts = read_data(filename)
times = np.array([t0+tstep*i for i in range(len(counts))])

def find_peaks(counts, times, spacing, cutoff): # Prints the locations of peaks in data
    
    calibration_factor = 0
    num = 0
    i = 0
    while i < len(counts):
        if counts[i] >= cutoff:
            calibration_factor = 0
            num = 0
            error = 0
            for j in range(spacing):
                if i+j < len(counts) and counts[i+j] >= cutoff:
                    calibration_factor += counts[i+j]*times[i+j]
                    error += tstep**2
                    num += counts[i+j]    
            print(f"{calibration_factor/num} +- {np.sqrt(error/num)}")
            i += spacing
        i += 1

find_peaks(counts, times, spacing, cutoff)