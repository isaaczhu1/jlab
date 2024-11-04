# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:00:12 2024

@author: evanb
"""

import matplotlib.pyplot as plt
import numpy as np
import math

def read_data(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    
    returns np array of all coincidences

    """
    
    
    f = open(f".data/{filename}.txt")
    for i in range(7):
        f.readline()
    
    fielded_data = []

    datum = f.readline()
    while not datum == "":
        data = datum.split()
        if data[10] == "1":
            fielded_data.append(np.array(data))
        datum = f.readline()
    f.close()
    return np.array(fielded_data)