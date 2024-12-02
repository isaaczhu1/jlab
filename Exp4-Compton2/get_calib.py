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
    'Cs137': 661.66
}

def gaussian(x, a, b, c):
    return a * np.exp(-(x - b)**2 / (2 * c**2))

def get_peak(counts, est_peak, plot=False, peak_name="", debug=False):
    '''
    Get the index of the peak in the counts array.

    Parameters:
        counts: the histogram array
        est_peak: an estimate of where the peak is, as an index in the counts array
        debug: whether to print debug information
    Returns:
        a dictionary with the following keys
            'mean': the mean of the peak as an index in the counts array
            'std': the standard deviation of the peak
            'amplitude': the amplitude of the peak
    '''
    # smooth out counts via convolution
    n_conv = 20
    counts = np.convolve(counts, np.ones(n_conv)/n_conv, mode='same')

    MAX_CNT = max(counts[est_peak - 30: est_peak + 30])

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
    peak_params, pcov = curve_fit(gaussian, peak_x, peak_counts, p0=[MAX_CNT, len(peak_counts) / 2, len(peak_counts)])
    peak_mean = int(peak_params[1] + min_peak_region_index)

    # print(f'Peak at {peak_mean}')

    if plot:
        plt.plot(min_peak_region_index + peak_x, gaussian(peak_x, *peak_params), color='red')
        plt.text(peak_mean*1.1, MAX_CNT, f'{peak_name} peak', fontsize=12, color='red')
        # plt.show()


    ret = {
        'mean': peak_mean, # the mean of the peak as an index in the counts array
        'mean_err': np.sqrt(np.diag(pcov))[1],
        'std': peak_params[2],
        'std_err': np.sqrt(np.diag(pcov))[2],
        'amplitude': peak_params[0],
        'amplitude_err': np.sqrt(np.diag(pcov))[0]
    }

    return ret


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

    plt.plot(counts, alpha=0.5)

    # Find the Na22 peak
    est_na22_peak = 500 + np.argmax(counts[500:])
    na22_data = get_peak(counts, est_na22_peak, plot=True, peak_name="Na22", debug=verbose)
    na22_mean = na22_data['mean']
    na22_std = na22_data['std']
    na22_amplitude = na22_data['amplitude']

    # plt.show()

    # Find the Cs137 peak
    est_cs137_peak = int(na22_mean * 1.15) + np.argmax(counts[int(na22_mean * 1.15):])
    cs137_data = get_peak(counts, est_cs137_peak, plot=True, peak_name="Cs137", debug=verbose)
    cs137_mean = cs137_data['mean']
    cs137_std = cs137_data['std']
    cs137_amplitude = cs137_data['amplitude']

    # plot the histogram with fitted peaks

    # get the name of the file. this consists of either "scatter" or "recoil" followed by the angle
    if "scatter" in filename:
        file_type = "scatter"
        angle = filename[7:-4]
    else:
        file_type = "recoil"
        angle = filename[6:-4]
    
    plt.title(f'Calibration of {file_type} at {angle}°', fontsize=16)
    plt.xlabel('MCA channel', fontsize=13)
    plt.ylabel('Counts', fontsize=13)
    plt.ylim(0, max(counts[100:]) * 1.1)
    plt.savefig(f'./images/calib_hists/{filename}.png')
    plt.clf()

    if verbose:
        print("Na22 peak energy:", na22_mean)
        print("Cs137 peak energy:", cs137_mean)

    # Calibrate the histogram.
    calib_slope = (known_energies['Cs137'] - known_energies['Na22']) / (cs137_mean - na22_mean)
    calib_intercept = known_energies['Na22'] - calib_slope * na22_mean

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
    angles = [30, 60, 90, 120, 135]
    marked_filenames = [f'scatter{angle}.Chn' for angle in angles] \
        + [f'recoil{angle}.Chn' for angle in angles]
    # get the calibrated x axes for each file, and store them in data/calib_x.pkl
    calib_x = {}
    for filename in marked_filenames:
        if "scatter" in filename:
            file_type = "scatter"
            angle = filename[7:-4]
        else:
            file_type = "recoil"
            angle = filename[6:-4]
        calib_x[filename] = get_calib_x(filename)
        # generate a plot of the calibrated histogram
        raw_counts = read_data(filename)['counts']
        plt.plot(calib_x[filename], raw_counts, alpha=0.5)
        plt.title(f'Calibrated {file_type} at {angle}°', fontsize=16)
        plt.xlabel('Energy (keV)', fontsize=13)
        plt.ylabel('Counts', fontsize=13)
        plt.savefig(f'./images/calibrated_hists/{filename}_calib.png')
        plt.clf()
        print(f'Calibrated {filename}')


    with open('./data/calib_x.pkl', 'wb') as f:
        pickle.dump(calib_x, f)
        print('Stored calibrated x axes in data/calib_x.pkl')