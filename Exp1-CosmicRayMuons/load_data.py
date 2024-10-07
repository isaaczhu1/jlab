'''

load data from binary file such as ./data/calibration.Chn

TODO: convert x-axis to time-difference

'''

import numpy as np
import matplotlib.pyplot as plt

def read_data(filename, left=16, right=256):
    with open(f'./data/{filename}.Chn', 'rb') as f:
        data = f.read()
        # the data type is INTEGERS
        counts = np.frombuffer(data, dtype=np.uint16)

    # print("counts.shape:",counts.shape)

    # must remove garbage stuff at the beginning and end
    counts = counts[left:-right]

    return counts

def calibrate_time(filename, calib_time):
    # return the averaged distance between peaks
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
    # print("diffs:", diffs)
    return calib_time / np.mean(diffs)

def remove_zeros_on_margin(counts):
    is_zero = True
    for i in range(len(counts)):
        if counts[i] != 0:
            break
    counts = counts[i:]
    return counts, i

if __name__ == '__main__':
    filename = 'lifetime_settings_calib'
    counts = read_data(filename)

    # print the indices where counts > 10000
    print(np.where(counts >= max(counts)))

    plt.plot(counts)
    plt.title(f'{filename}')
    plt.xlabel('Voltage bin')
    plt.ylabel('Counts')
    plt.savefig(f'./images/{filename}.png')