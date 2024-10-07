<<<<<<< HEAD
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
=======
from load_data import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = read_data('lifetime_weekend')
data, num_leading_zeros = remove_zeros_on_margin(data)

# bin the data
BIN_SIZE = 20
CALIB_TIMESCALE = 1.28 # in mus

TIME_PER_BIN = BIN_SIZE * calibrate_time('lifetime_settings_calib', CALIB_TIMESCALE) # in mus
print(f"TIME_PER_BIN: {TIME_PER_BIN:.5f} mus")
data = BIN_SIZE*np.mean(data[:len(data)//BIN_SIZE * BIN_SIZE].reshape(-1, BIN_SIZE), axis=1)
print("num bins after binning:", data.shape)


# fit an exponential decay curve
# y = A * exp(-x/tau) + C
# where A is the amplitude, tau is the lifetime, and C is the offset

def exponential_decay(x, A, tau, C):
    return A * np.exp(-x / tau) + C

def pure_exponential_decay(x, A, tau):
    return A * np.exp(-x / tau)

def fit_exponential_decay(x, y):
    popt, pcov = curve_fit(exponential_decay, x, y, p0=[max(y), 2/TIME_PER_BIN, 0])
    return popt, pcov

def fit_pure_exponential_decay(x, y):
    popt, pcov = curve_fit(pure_exponential_decay, x, y, p0=[max(y), 2/TIME_PER_BIN])
    return popt, pcov

x = np.arange(len(data))
y = data

USE_BACKGROUND = False

if USE_BACKGROUND:
    popt, pcov = fit_exponential_decay(x, y)
    # compute the chi-squared
    y_pred = exponential_decay(x, *popt)
    print((y - y_pred)**2 / y_pred)
    chi_squared = np.sum((y - y_pred)**2 / y_pred)
    dof = len(y) - len(popt)
else:
    popt, pcov = fit_pure_exponential_decay(x, y)
    # compute the chi-squared
    y_pred = pure_exponential_decay(x, *popt)
    chi_squared = np.sum((y - y_pred)**2 / y_pred)
    dof = len(y) - len(popt)

print("popt",popt)
# the conversion factor is (BIN_SIZE * 10)/4096 ms per bin
tau = popt[1] * TIME_PER_BIN
tau_err = np.sqrt(np.diag(pcov))[1] * TIME_PER_BIN

x_times = x * TIME_PER_BIN

SHIFT_TIMES = True

if SHIFT_TIMES:
    x_times = x_times + (num_leading_zeros // BIN_SIZE) * TIME_PER_BIN

print(f'Lifetime: {tau:.3f} +/- {tau_err:.3f} $\\mu$s')

# plt.scatter(x, y, label='data', alpha=0.5)
plt.scatter(x_times, y, label='data', alpha=0.5)
# plt.plot(x, y_pred, label='fit', color='orange')
plt.plot(x_times, y_pred, label='exponential fit', color='#ff7f0e')
plt.title('Lifetime of muons')
plt.xlabel('Decay Time ($\\mu$s)')
plt.ylabel('Binned Counts')

plt.text(0.5, 0.6, f'$\\chi^2 / \\text{{dof}} = {chi_squared:.1f} / {dof}$', transform=plt.gca().transAxes)
plt.text(0.5, 0.5, f'$\\tau = {tau:.3f} \pm {tau_err:.3f}$ $\\mu$s', transform=plt.gca().transAxes)

plt.legend()
plt.savefig('./images/lifetime_weekend.png')
>>>>>>> f7f4df3b3d7a4ac1225840c1bdb7e15f4202c0e1
