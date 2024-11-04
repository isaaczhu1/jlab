'''

an example of how to load data from data/nursery.h5

'''

import h5py

with h5py.File('data/nursery.h5', 'r') as f:
    group_name = 'data_nursery_TRPM_FA2024_trpm-the_dudes_'
    for key in f[group_name].keys():
        print(f"key has shape:")
        print(f[group_name][key][:].shape)