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

def get_tof(filename, generate_plot=False):
    # filename = "150cmlongrun"
    bin_width = 1/79.15 # in ns TODO verify?

    data = read_data(filename)
    centers = [i*bin_width for i in range(len(data))]

    BIN_WIDTH = 32

    valid_data = np.zeros(math.floor(len(data)/BIN_WIDTH))
    valid_centers = np.zeros(math.floor(len(data)/BIN_WIDTH))

    for i in range(len(data)):
        valid_data[math.floor(i/BIN_WIDTH)] += data[i]
        valid_centers[math.floor(i/BIN_WIDTH)] += centers[i]/BIN_WIDTH

    # remove any bin where the count is 0. the result is called valid_data
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
    # normal = lambda x, m, s, A, C: A*np.exp(-(x-m)**2/(2*s**2))+C
    # lorentz = lambda x, m, s, A, C: A/((x-m)**2+s**2)+C

    # initial_guess = [40, 10, 200, 0.1]

    est_mean = valid_centers[np.argmax(valid_data)]
    est_std = 10
    est_A = valid_data[np.argmax(valid_data)]
    est_alpha = 0.1
    initial_guess = [est_mean, est_std, est_A, est_alpha, 0]

    # we will fit a SKEWED NORMAL distribution to (valid_centers, valid_data)
    # the skewed normal distribution is defined as:
    # f(x) = 2*phi(x)*Phi(alpha*x)
    # where phi(x) is the standard normal distribution and Phi(x) is the cumulative distribution function of the standard normal distribution
    # Phi(x) = 1/2*(1+erf(x/sqrt(2)))

    def skew_normal(x, m, s, A, alpha, C):
        return A*2*1/np.sqrt(2*np.pi)*np.exp(-1/2*(x-m)**2/s**2)*1/2*(1+scipy.special.erf(alpha*(x-m)/np.sqrt(2)))+C

    # initial_guess = [40, 10, 200, 0.1, 0]

    [mean, std, amp, alpha, c], cov = scipy.optimize.curve_fit(skew_normal, valid_centers, valid_data, p0=initial_guess,sigma=np.sqrt(valid_data))
    # print("params", [mean, std, amp, alpha, c])

    # error = np.sqrt(cov[0][0])

    y_pred = np.array([skew_normal(x,mean, std, amp, alpha, c) for x in valid_centers])

    # compute the chi squared value
    chi = sum([(valid_data[i]-y_pred[i])**2/y_pred[i] for i in range(len(y_pred))])
    dof = len(valid_data) - 5

    # we're going to define the flight time as the PEAK OF THE SKEWED NORMAL DISTRIBUTION
    # find the index of the maxmium value of y_pred
    fit_peak = valid_centers[np.argmax(y_pred)]
    peak_error = np.sqrt(cov[0][0])
    if generate_plot:
        # plot the data and the fit
        plt.scatter(valid_centers, valid_data,alpha=0.5,marker=".")
        plt.plot(valid_centers, y_pred)
        # plt.annotate(f"Mean: {np.round(mean,2)}$\pm${np.round(error,2)} ns", (10,25))
        plt.annotate(f"Peak: {np.round(fit_peak,2)}$\pm${np.round(peak_error,2)} ns", (7, est_A/2))
        plt.annotate(f"$\chi^2 / dof$: {np.round(chi)} / {dof}", (7, est_A/3))
        plt.title(f'{filename}')
        plt.xlabel('Time (ns)')
        plt.ylabel('Counts')
        plt.savefig(f'./images/{filename}_tof_fit.png')

    return fit_peak, peak_error


if __name__ == '__main__':
    filename = '290cm_final'
    print(get_tof(filename, generate_plot=True))