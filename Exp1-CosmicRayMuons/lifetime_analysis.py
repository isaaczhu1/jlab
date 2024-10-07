from load_data import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = read_data('lifetime_weekend')

print("original num bins:", data.shape)
# plt.plot(data)
# plt.show()

# assert False

data, num_leading_zeros = remove_zeros_on_margin(data)
print("num bins after removing dead zone with zeros:", data.shape)

# bin the data
BIN_SIZE = 20
TIME_PER_BIN = BIN_SIZE * 10 / 4096 # in mus
data = 20*np.mean(data[:len(data)//BIN_SIZE * BIN_SIZE].reshape(-1, BIN_SIZE), axis=1)
print("num bins after binning:", data.shape)


# fit an exponential decay curve
# y = A * exp(-x/tau) + C
# where A is the amplitude, tau is the lifetime, and C is the offset

def exponential_decay(x, A, tau, C):
    return A * np.exp(-x / tau) + C

def pure_exponential_decay(x, A, tau):
    return A * np.exp(-x / tau)

def fit_exponential_decay(x, y):
    popt, pcov = curve_fit(exponential_decay, x, y, p0=[max(y), 100, 0])
    return popt, pcov

def fit_pure_exponential_decay(x, y):
    popt, pcov = curve_fit(pure_exponential_decay, x, y, p0=[max(y), 100])
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
    x_times = x_times + num_leading_zeros // 20 * TIME_PER_BIN

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