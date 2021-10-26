#!/usr/bin/python3

import os
from collections import OrderedDict
import numpy as np
import h5py

def add_to_key_dict(key_dict, tokens):
    fname       = tokens[0]
    key_hdf5    = tokens[1]
    key_import  = tokens[2]
    description = tokens[3]
    if fname in key_dict.keys(): pass
    else: key_dict[fname] = []
    key_dict[fname].append((key_hdf5, key_import, description))

def main():
    # Get list of tallies
    with open('tally_list.txt', 'r') as r:
        tally_list = [int(x) for x in r.read().split()]

    # Get list of expected datasets
    key_dict = OrderedDict()
    with open('../datasets.txt', 'r') as r: lines = r.readlines()
    for line in lines:
        if line.startswith('+'): continue
        tokens = [x.strip() for x in line.split('|')[1:-1]]
        if not tokens[0].endswith('h5'): continue
        if '{T}' in line:
            for tally in tally_list:
                line_new = line.replace('{T}', str(tally))
                tokens_new = [x.strip() for x in line_new.split('|')[1:-1]]
                add_to_key_dict(key_dict, tokens_new)
        else: add_to_key_dict(key_dict, tokens)

    # Print header
    hriz_line = ('+%s+%s+%s+%s+%s+%s+%s+%s+' %
                 ('-' * 29, '-' * 27, '-' * 26, '-' * 24, '-' * 10, '-' * 13,
                  '-' * 13, '-' * 40))
    print(hriz_line)
    print('| %-27s | %-25s | %-24s | %-22s | %-8s | %11s | %11s | %-38s |' %
          ('Filename', 'HDF5 key', 'Imported key', 'Shape', 'Type',
           'Num Nonzero', 'Num Zero', 'Description'))
    print(hriz_line)

    # Print information about each dataset
    for fname in key_dict.keys():
        if not os.path.exists(fname): continue
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
            print(('| %-27s | %-25s | %-24s | %-22s | %-8s | %11s | %11s '
                   '| %-38s |') %
                  (fname, key_hdf5, key_import, val.shape, dtype_str,
                   num_nonzero, num_zero, description))
        print(hriz_line)

if __name__ == '__main__': main()
