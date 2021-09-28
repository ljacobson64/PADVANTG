#!/usr/bin/python3

import os
from collections import OrderedDict
import numpy as np
import h5py

def load_hdf5_data(group):
    data = {}
    for key in group:
        item = group.get(key)
        if isinstance(item, h5py.Dataset): data[key] = np.array(item[()])
        else:                              data[key] = load_hdf5_data(item)
    return data

def main():
    key_dict = OrderedDict()
    with open('datasets.txt', 'r') as r: lines = r.readlines()
    for line in lines:
        if line.startswith('+'): continue
        tokens = [x.strip() for x in line.split('|')[1:-1]]
        if not tokens[0].endswith('h5'): continue
        fname       = tokens[0]
        key_hdf5    = tokens[1]
        key_import  = tokens[2]
        description = tokens[3]
        if fname in key_dict.keys(): pass
        else: key_dict[fname] = []
        key_dict[fname].append((key_hdf5, key_import, description))

    hriz_line = ('+%s+%s+%s+%s+%s+%s+%s+%s+' %
                 ('-' * 36, '-' * 27, '-' * 20, '-' * 23, '-' * 10, '-' * 13,
                  '-' * 13, '-' * 40))
    print(hriz_line)
    print('| %-34s | %-25s | %-18s | %-21s | %-8s | %11s | %11s | %-38s |' %
          ('Filename', 'HDF5 key', 'Imported key', 'Shape', 'Type',
           'Num Nonzero', 'Num Zero', 'Description'))
    print(hriz_line)

    for fname in key_dict.keys():
        hf = h5py.File(fname, 'r')
        for key_hdf5, key_import, description in key_dict[fname]:
            val = hf[key_hdf5][()]
            dtype_str = val.dtype
            if key_hdf5 == 'mixtable': dtype_str = 'i4,i4,f8'
            if key_hdf5 in ['mixtable', 'mat_names']:
                num_nonzero = val.shape[0]
                num_zero = 0
            else:
                num_nonzero = np.count_nonzero(val)
                num_zero = np.count_nonzero(val == 0)
            print(('| %-34s | %-25s | %-18s | %-21s | %-8s | %11s | %11s '
                   '| %-38s |') %
                  (fname, key_hdf5, key_import, val.shape, dtype_str,
                   num_nonzero, num_zero, description))
        print(hriz_line)

if __name__ == '__main__': main()
