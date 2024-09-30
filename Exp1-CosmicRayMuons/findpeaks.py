"""
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

v0 = 0 #Starting bin voltage
vstep = 10/len(counts)  # Bin size
cutoff = 20 # Number of counts required to be considered in max

spacing = 100 # Must be under the spacing between peaks in volts.

counts = read_data(filename)
voltages = np.array([v0+vstep*i for i in range(len(counts))])

def find_peaks(counts, voltages, spacing, cutoff): # Prints the locations of peaks in data
    
    calibration_factor = 0
    error = 0
    num = 0
    i = 0
    while i < len(counts):
        if counts[i] >= cutoff:
            calibration_factor = 0
            num = 0
            for j in range(spacing):
                if i+j < len(counts) and counts[i+j] >= cutoff:
                    calibration_factor += counts[i+j]*voltages[i+j]
                    num += counts[i+j]    
            print(calibration_factor/num)
            i += spacing
        i += 1
