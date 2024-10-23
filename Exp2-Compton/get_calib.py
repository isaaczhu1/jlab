'''

There are various peaks in the calibration files 
['scatterNaCalib',
'recoilNaCalib',
'scatterBaCalib',
'recoilBaCalib',
'scatterCsCalib',
'recoilCsCalib',]

(scatter and recoil are two separate detectors. We used Na, Ba, and Cs as calibration sources for each detector.)           

There are known energies for each of the peaks in the calibration files. 
We will use this to compute a linear calibration function (from channel number to energy) for each detector.

'''

from load_data import *

known_energies = {
    'Na' : [511,],
    'Ba' : [81, 302.85, 356.02,],
    'Cs' : [661.66,],
}

# The calibration function will be of the form Energy (in keV) = factor * Channel number + intercept
# As a first approximation, we estimate factor ~ 0.17 and intercept ~ 0
naive_factor = 0.17
naive_intercept = 0.0

def get_calib(filename, known_energies, factor, intercept):
    '''
    Compute a pair of the form (channel number, known energy) for each peak in the calibration file.
    Use this to compute the calibration factor and intercept.
    '''
    raise NotImplementedError
    counts = read_data(filename)

    # compute the estimated channel number for each peak, using the naive factor and intercept
    peaks = []
    for energy in known_energies[filename[7:9]]:
        peaks.append(energy / factor - intercept)