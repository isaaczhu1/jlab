'''

For each file e.g. scatter100.Chn,
there's another file scatter100marked.Chn
which is also a histogram, but with large Na22 and Cs137 peaks on them.
These peaks correspond to energies of 511 keV and 662 keV respectively.
We can use these peaks to calibrate the energy of the original histogram in scatter100.Chn.

'''

from load_data import *
from scipy.optimize import curve_fit

known_energies = {
    'Na22': 511,
    'Cs137': 662
}

def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (2 * c**2))

def get_peak(counts, est_peak, debug=False):
    '''
    Given a file and an estimated peak, returns the index of the peak by fitting a gaussian distribution around the peak.
    '''
    # smooth out counts via convolution
    orig_counts = counts
    n_conv = 20
    counts = np.convolve(counts, np.ones(n_conv)/n_conv, mode='same')
    if debug:
        plt.clf()
        plt.plot(counts, alpha=0.5, label='smoothed')
        plt.plot(orig_counts, alpha=0.5, label='original')
        plt.legend()
        plt.show()
        plt.clf()
    MAX_CNT = max(counts[est_peak - 50: est_peak + 50])

    # Get a region around the peak
    def is_peak(i):
        return counts[i] > MAX_CNT // 2
    peak_region = [0 for _ in range(len(counts))]
    peak_region[est_peak] = 1

    i = est_peak
    while i > 0 and is_peak(i):
        i -= 1
        peak_region[i] = 1
    i = est_peak
    while i < len(counts) - 1 and is_peak(i):
        i += 1
        peak_region[i] = 1

    min_peak_region_index = min([i for i in range(len(peak_region)) if peak_region[i] == 1])

    peak_counts = counts[np.array(peak_region) == 1]
    peak_x = np.arange(len(peak_counts))
    peak_params, _ = curve_fit(gaussian, peak_x, peak_counts, p0=[MAX_CNT, len(peak_counts) / 2, len(peak_counts)])
    peak_mean = int(peak_params[1] + min_peak_region_index)


    if debug:
        print(f"Estimated peak index: {est_peak}, actual peak index: {peak_mean}")
        plt.clf()
        plt.plot(counts)
        plt.plot(peak_x + min_peak_region_index, gaussian(peak_x, *peak_params))
        # plt.plot(peak_x, peak_counts)
        # plt.plot(peak_x, gaussian(peak_x, *peak_params))
        plt.title("Peak fit around index {}".format(est_peak))
        plt.show()

    return int(peak_mean)

def get_calib_x(filename, verbose=False):
    '''
    Calibrates based off Na22 and Cs137 peaks.
    Parameters:
        file: the name of the file to read from, which encodes a histogram
    Returns:
        the calibrated histogram array.
    '''
    calib_file = filename[:-4] + 'marked.Chn'
    counts = read_data(calib_file)['counts']  # a numpy array

    # Find the Na22 peak
    est_na22_peak = 1000 + np.argmax(counts[1000:])
    na22_mean = get_peak(counts, est_na22_peak, debug=verbose)

    # Find the Cs137 peak
    est_cs137_peak = int(na22_mean * 1.15) + np.argmax(counts[int(na22_mean * 1.15):])
    cs137_mean = get_peak(counts, est_cs137_peak, debug=verbose)

    if verbose:
        print("Na22 peak energy:", na22_mean)
        print("Cs137 peak energy:", cs137_mean)

    # Calibrate the histogram.
    calib_slope = (known_energies['Cs137'] - known_energies['Na22']) / (cs137_mean - na22_mean)
    calib_intercept = known_energies['Na22'] - calib_slope * na22_mean

    if verbose:
        print("Calibration: y = {}x + {}".format(calib_slope, calib_intercept))

    # Now we get the calibrated x axis values.
    calib_counts = read_data(filename)['counts']
    calib_x = np.arange(len(calib_counts))
    calib_x = calib_slope * calib_x + calib_intercept
    return calib_x

def inv_calib_x(orig_x, calib_x, val):
    '''
    orig_x and calib_x are the original and calibrated lists of x values respectively.
    val is the value to calibrate.
    '''
    slope = (calib_x[-1] - calib_x[0]) / (orig_x[-1] - orig_x[0])
    intercept = calib_x[0] - slope * orig_x[0]
    return int((val - intercept) / slope)

if __name__ == '__main__':
    calibrated_x = get_calib_x('scatter100.Chn', verbose=False)
    plt.plot(calibrated_x, read_data('scatter100marked.Chn')['counts'])
    plt.show()