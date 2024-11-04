# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:28:06 2024

@author: evanb
"""
from readfile import read_data
import numpy as np
import matplotlib.pyplot as plt

# Extract file name
filename = "" 

data = read_data(filename)

start_time = 0 # in seconds after midnight

def extract_times(clock):
    """
    

    Parameters
    ----------
    clock : string
        the clock time of the data.

    Returns
    -------
    integer
        seconds after midnight that the detection was made.

    """
    time_splits = clock.split(":")
    time = (3600*int(time_splits[0])+60*int(time_splits[1])+int(time_splits[2]))+start_time
    return time % (24*3600)


# Get times after midnight for every single coincidence
times = [] 

for datum in data:
    times.append(extract_times(datum[1]))

plt.hist(times) # See if there is some correlation in time
    


