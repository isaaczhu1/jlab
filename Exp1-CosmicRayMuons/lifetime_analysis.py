'''
author: isaaczhu
'''

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
# print("num bins after binning:", data.shape)


# fit an exponential decay curve y = A * exp(-x/tau) + C

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

USE_BACKGROUND = False # whether to fit an exponential with a background term

if USE_BACKGROUND:
    popt, pcov = fit_exponential_decay(x, y)
    # compute the chi-squared
    y_pred = exponential_decay(x, *popt)
    # print((y - y_pred)**2 / y_pred)
    chi_squared = np.sum((y - y_pred)**2 / y_pred)
    dof = len(y) - len(popt)
else:
    popt, pcov = fit_pure_exponential_decay(x, y)
    # compute the chi-squared
    y_pred = pure_exponential_decay(x, *popt)
    chi_squared = np.sum((y - y_pred)**2 / y_pred)
    dof = len(y) - len(popt)

# print(popt)

# compute tau (the mean lifetime) in microseconds
tau = popt[1] * TIME_PER_BIN
tau_err = np.sqrt(np.diag(pcov))[1] * TIME_PER_BIN

x_times = x * TIME_PER_BIN # for plotting

SHIFT_TIMES = True
if SHIFT_TIMES:
    # account for the fact that the data is shifted by the number of leading zeros
    x_times = x_times + (num_leading_zeros // BIN_SIZE) * TIME_PER_BIN

print(f'Lifetime: {tau:.3f} +/- {tau_err:.3f} $\\mu$s')

plt.scatter(x_times, y, label='data', alpha=0.5)
plt.plot(x_times, y_pred, label='exponential fit', color='#ff7f0e')
plt.text(0.5, 0.6, f'$\\chi^2 / \\text{{dof}} = {chi_squared:.1f} / {dof}$', transform=plt.gca().transAxes)
plt.text(0.5, 0.5, f'$\\tau = {tau:.3f} \pm {tau_err:.3f}$ $\\mu$s', transform=plt.gca().transAxes)


plt.title('Lifetime of muons', fontsize=15)
plt.xlabel('Decay Time ($\\mu$s)', fontsize=13)
plt.ylabel('Binned Counts', fontsize=13)
plt.legend()
plt.savefig('./images/lifetime_weekend.png')