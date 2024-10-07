# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 09:36:25 2024

@author: evanb
"""

from load_data import read_data
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib
import math

filename = "bad150cm"
bin_width = 1/79.15 # in ns

data = read_data(filename)
centers = [i*bin_width for i in range(4096)]

BIN_WIDTH = 32

valid_data = np.zeros(math.floor(len(data)/BIN_WIDTH))
valid_centers = np.zeros(math.floor(len(data)/BIN_WIDTH))

for i in range(len(data)):
    valid_data[math.floor(i/BIN_WIDTH)] += data[i]
    valid_centers[math.floor(i/BIN_WIDTH)] += centers[i]/BIN_WIDTH


bad_index = []
for i in range(len(valid_data)):
    if valid_data[i] == 0:
        bad_index.append(i)
for i in range(len(bad_index)-1, -1, -1):
    valid_data = np.delete(valid_data,bad_index[i])
    valid_centers = np.delete(valid_centers,bad_index[i])

valid_data = np.delete(valid_data,len(valid_data)-1)
valid_centers = np.delete(valid_centers,len(valid_data)-1)
valid_data = np.array(valid_data)
valid_centers = np.array(valid_centers)
normal = lambda x, m, s, A, C: A*np.exp(-(x-m)**2/(2*s**2))+C
lorentz = lambda x, m, s, A, C: A/((x-m)**2+s**2)+C



initial_guess = [40, 10, 200,0]
means = []


[mean, std, amp,c], cov = scipy.optimize.curve_fit(normal, valid_centers, valid_data, p0=initial_guess,sigma=np.sqrt(valid_data))

error = np.sqrt(cov[0][0])

lor = np.array([normal(x,mean, std, amp,c) for x in valid_centers])

chi = sum([(valid_data[i]-lor[i])**2/lor[i] for i in range(len(lor))])


plt.scatter(valid_centers, valid_data,alpha=0.5,marker=".")
plt.plot(valid_centers, lor)
plt.annotate(f"Mean: {np.round(mean,2)}$\pm${np.round(error,2)} ns\n$\chi^2$: {np.round(chi)} / {len(valid_data - 3)}", (10,25))
plt.title(f'{filename}')
plt.xlabel('Time (ns)')
plt.ylabel('Counts')
