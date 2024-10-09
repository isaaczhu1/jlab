# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 09:20:02 2024

@author: evanb
"""
import numpy as np
import random

# Supply the the data and their standard deviations

def lin_regression_with_error(x, y, xerr, yerr, trials):
    y_values = np.array(y)
    x_values = np.array(x)
    y_err = np.array(xerr)
    x_err = np.array(yerr)

    trials = trials # How many monti carlo runs to make

    # Calculate best fit line
    m, b = np.polyfit(x_values, y_values,1)

    # Run monti carlo on data
    m_trials = []
    b_trials = []
    for i in range(trials):
        y_trial = np.array(y_values)
        x_trial = np.array(x_values)
        for i in range(len(x_values)):
            x_trial[i] += x_err[i]*random.gauss(0,1)
            y_trial[i] += y_err[i]*random.gauss(0,1)
        mt, bt = np.polyfit(x_trial, y_trial,1)
        m_trials.append(mt)
        b_trials.append(bt)
    m_trials = np.array(m_trials)
    b_trials = np.array(b_trials)
    print(len(b_trials))
    # Use data to give standard deviation
    m_std = np.std(m_trials)
    b_std = np.std(b_trials)
    return m, b, m_std, b_std
    
    