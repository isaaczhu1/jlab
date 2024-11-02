import struct
import numpy as np
import matplotlib.pyplot as plt
import json


def read_data(filename, left=16, right=256):
    '''
    Reads data from filename.Chn (which encodes a histogram)

    Returns:
        a dictionary with the following keys
        - counts: the counts in each bin
        - time: the time the data was taken for
    '''
    # Chn constants
    HEADER_SIZE = 32
    CHANNELS_PER_RECORD = 8
    with open(f'./data/{filename}', 'rb') as f:
        # Read the first 32 bytes of header data
        header_data = f.read(HEADER_SIZE)
        # Unpack the header data (Fortran 'INTEGER*2' is 'int16' in Python and 'INTEGER*4' is 'int32')
        TYPE, MCA, SEG, SRTSEC, RLTIME, LVETME, SRTDTE, SRTTME, STRTCH, LNGTDT = struct.unpack('h h h 2s i i 8s 4s h h', header_data)
        data = f.read()
        # the data type is INTEGERS
        counts = np.frombuffer(data, dtype=np.uint16)

    # must remove garbage stuff at the beginning and end
    counts = counts[left:-right]

    # make sure counts is writable
    counts = counts.copy()

    ret = {
        'counts': counts,
        'time': LVETME / 50,
    }

    return ret

def calibrate_time(filename, calib_time, verbose=False):
    '''
    Parameters:
        filename: the name of the file to read from, which encodes a histogram
        calib_time: the distance between peaks in the calibration file
    Returns:
        the time per bin
    '''
    counts = read_data(filename)['counts']
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

def gen_plot(filename, binsize=1, scale_from=0):
    counts = read_data(filename)['counts']
    # find the largest value with index greater than scale_from
    max_val = max(counts[scale_from:])
    rebined = rebin(counts, binsize)

    plt.plot(rebined)
    plt.title(f'{filename}')
    plt.xlabel('Voltage bin', fontsize=13)
    plt.ylabel('Counts', fontsize=13)
    plt.ylim(0, max_val*1.5)
    plt.savefig(f'./images/{filename}.png')
    plt.clf()

if __name__ == '__main__':
    # generate all the raw histograms
    angles = [30, 60, 90, 120, 140]
    filenames = [f'goode/scatter{angle}.Chn' for angle in angles] \
        + [f'goode/scatter{angle}marked.Chn' for angle in angles] \
        + [f'goode/recoil{angle}.Chn' for angle in angles] \
        + [f'goode/recoil{angle}marked.Chn' for angle in angles]
    for filename in filenames:
        gen_plot(filename, binsize=1, scale_from=100)
        print(f'Generated plot for {filename}')

    # store the live times of all the filenames in a json
    live_times = {filename: read_data(filename)['time'] for filename in filenames}
    with open('./data/live_times.json', 'w') as f:
        json.dump(live_times, f)
    print('Stored live times in live_times.json')