'''

Load the data, generate plots of raw histograms, and store the data in a pickle file

'''


import struct
import numpy as np
import matplotlib.pyplot as plt
import pickle
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

    # rebin counts by a factor of 2
    counts = np.array([sum(counts[i:i+2]) for i in range(0, len(counts), 2)])

    # make sure counts is writable
    counts = counts.copy()

    ret = {
        'counts': counts,
        'time': LVETME / 50,
    }

    return ret


def gen_plot(filename, binsize=1, scale_from=0):
    counts = read_data(filename)['counts']
    # find the largest value with index greater than scale_from
    max_val = max(counts[scale_from:])

    if "scatter" in filename:
        file_type = "scatter"
    else:
        file_type = "recoil"
    if "marked" in filename:
        marked = True
    else:
        marked = False
    # angle is the digits in the file name
    angle = ""
    for char in filename:
        if char.isdigit():
            angle += char

    title = f'Raw data for {file_type} at {angle}°'
    if marked:
        title += " (marked)"

    plt.plot(counts)
    plt.title(title, fontsize=16)
    plt.xlabel('MCA channel', fontsize=13)
    plt.ylabel('Counts', fontsize=13)
    plt.ylim(0, max_val*1.5)
    # plt.show()
    plt.savefig(f'./images/raw_hists/{filename}.png')
    plt.clf()

if __name__ == '__main__':
    # generate all the raw histograms
    angles = [30, 60, 90, 120, 135]
    filenames = [f'scatter{angle}.Chn' for angle in angles] \
        + [f'scatter{angle}marked.Chn' for angle in angles] \
        + [f'recoil{angle}.Chn' for angle in angles] \
        + [f'recoil{angle}marked.Chn' for angle in angles]
    for filename in filenames:
        gen_plot(filename, binsize=1, scale_from=100)
        print(f'Generated plot for {filename}')

    # store the following information in a pickle file: the the live times for each file, and a numpy array of the counts for each file
    data = {}
    for filename in filenames:
        data[filename] = read_data(filename)
    # with open('./data/data.pickle', 'wb') as f:
    #     pickle.dump(data, f)
    print('Stored data in data.pickle')

    # also store the live times in a json file
    live_times = {filename: data[filename]['time'] for filename in filenames}
    with open('./data/live_times.json', 'w') as f:
        json.dump(live_times, f, indent=4)

    # store the peak information in a json file
    peak_info = {}
    for filename in filenames:
        peak_info[filename] = {
            'peak': None,
            'num_counts': None,
        }
    with open('./data/peak_info.json', 'w') as f:
        json.dump(peak_info, f, indent=4)

    print('Stored live times in live_times.json')