from load_data import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = read_data('lifetime_weekend')

print(data.shape)
# plt.plot(data)
# plt.show()

# assert False

data = remove_zeros_on_margin(data)
print(data.shape)

# bin the data
bin_size = 20
data = np.mean(data[:len(data)//bin_size * bin_size].reshape(-1, bin_size), axis=1)
print(data.shape)


# fit an exponential decay curve
# y = A * exp(-x/tau) + C
# where A is the amplitude, tau is the lifetime, and C is the offset

def exponential_decay(x, A, tau, C):
    return A * np.exp(-x / tau) + C

x = np.arange(len(data))
y = data

popt, pcov = curve_fit(exponential_decay, x, y, p0=[max(y), 100, min(y)])
# print(popt)

# the conversion factor is (bin_size * 10)/4096 ms per bin
tau= popt[1] * (bin_size * 10) / 4096
tau_err = np.sqrt(np.diag(pcov))[1] * (bin_size * 10) / 4096
print(f'Lifetime: {tau:.3f} +/- {tau_err:.3f} ms')

plt.plot(x, y, label='data')
plt.plot(x, exponential_decay(x, *popt), label='fit')
plt.title('Lifetime of muons')
plt.xlabel('Time (ms)')
plt.ylabel('Binned Counts')

plt.text(0.5, 0.5, f'$\\tau = {tau:.3f} \pm {tau_err:.3f}$ ms', transform=plt.gca().transAxes)

plt.legend()
plt.savefig('./images/lifetime_weekend.png')