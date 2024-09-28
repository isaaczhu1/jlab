'''

load data from binary file such as ./data/calibration.Chn

TODO: convert x-axis to time-difference

'''

import numpy as np
import matplotlib.pyplot as plt

def read_data(filename):
    with open(f'./data/{filename}.Chn', 'rb') as f:
        data = f.read()
        # the data type is INTEGERS
        counts = np.frombuffer(data, dtype=np.uint16)

    print("counts.shape:",counts.shape)

    # must remove garbage stuff at the beginning and end
    counts = counts[16:-256]

    return counts

filename = 'test_run'
counts = read_data(filename)


plt.plot(counts)
plt.title(f'{filename}')
plt.xlabel('Voltage bin')
plt.ylabel('Counts')
plt.savefig(f'./images/{filename}.png')