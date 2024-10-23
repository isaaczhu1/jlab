import numpy as np
import matplotlib.pyplot as plt

def read_data(filename, left=16, right=256):
    '''
    Reads data from filename.Chn (which encodes a histogram) and returns the counts as a 1D numpy array
    '''
    with open(f'./data/{filename}.Chn', 'rb') as f:
        data = f.read()
        # the data type is INTEGERS
        counts = np.frombuffer(data, dtype=np.uint16)

    # must remove garbage stuff at the beginning and end
    counts = counts[left:-right]

    return counts

def calibrate_time(filename, calib_time, verbose=False):
    '''
    Parameters:
        filename: the name of the file to read from, which encodes a histogram
        calib_time: the distance between peaks in the calibration file
    Returns:
        the time per bin
    '''
    counts = read_data(filename)
    MAX_CNT = max(counts)
    def is_peak(i):
        return counts[i] > MAX_CNT // 4
    peaks = []
    for i in range(len(counts)):
        if is_peak(i):
            if len(peaks) == 0 or i - peaks[-1] > 10:
                peaks.append(i)
    diffs = np.diff(peaks)
    if verbose:
        print("number of bins between peaks:", diffs)
    return calib_time / np.mean(diffs)

def remove_zeros_on_margin(counts):
    '''
    For some of the data (in particular, the muon lifetime data), 
    there are zeros on the left margin of the histogram because the apparatus
    cannot pickup sufficiently short-lived muon decays.

    This function removes the zeros on the left margin of the histogram (truncates the array of counts).
    Also returns the number of zeros removed.
    '''
    for i in range(len(counts)):
        if counts[i] != 0:
            break
    counts = counts[i:]
    return counts, i

def rebin(counts, binsize):
    '''
    

    Parameters
    ----------
    counts: Data
    binsize : Number of data points placed into same bin. Accomplished by adding counts

    Returns
    Rebinned data, clips final incomplete bin

    '''
    rebin = np.zeros(len(counts)//binsize)
    for i in range(len(rebin)):
        for j in range(binsize):
            rebin[i] += counts[i*binsize+j]
    return rebin

if __name__ == '__main__':
    filename = 'scatter100'
    counts = read_data(filename)
    rebined = rebin(counts, 8)

    # print the indices where counts > 10000
    # print(np.where(counts >= max(counts)))

    plt.plot(rebined)
    plt.title(f'{filename}')
    #plt.title('Scattered Data', fontsize=15)
    plt.xlabel('Voltage bin', fontsize=13)
    plt.ylabel('Counts', fontsize=13)
    plt.savefig(f'./images/{filename}.png')
