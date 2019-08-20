#!/usr/bin/python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import h5py

def main():
    data = read_hdf5()
    print('')
    print('%-18s  %-7s  %-21s  %9s' % ('Key', 'Type', 'Shape', 'Size'))
    #sorted_keys = sorted(data)
    sorted_keys = sorted(sorted(data), key=lambda x: data[x].size, reverse=True)
    for key in sorted_keys:
        print('%-18s  %-7s  %-21s  %9s' %
              (key, data[key].dtype, data[key].shape, data[key].size))

def write_pickle(data, key):
    fname = 'pickles/%s.npy' % (key)
    print('Writing %s to file' % (fname))
    np.save(fname, data[key])

def read_hdf5():
    # Read HDF5 file containing forward flux
    fname = 'fwd_solution/denovo_forward.h5'
    print('Reading %s from file' % (fname))
    hf = h5py.File(fname,'r')

    # Load forward data
    data = {}
    for key in hf:
        if key in ['angles', 'quadrature_weights'] or key.startswith('mesh'):
            new_key = str(key)
        elif key in ['angular_flux']: new_key = 'angular_flux_fwd'
        else: continue
        data[new_key] = hf.get(key)[()]
        if data[new_key].dtype == np.int32:
            data[new_key] = np.array(data[new_key], dtype=np.int64)

    # Read HDF5 file containing adjoint flux
    fname = 'adj_solution/denovo_adjoint.h5'
    print('Reading %s from file' % (fname))
    hf = h5py.File(fname, 'r')

    # Load adjoint data
    data['angular_flux_adj'] = hf.get('angular_flux')[()]

    # Move energy axis to end
    for key in ['angular_flux_fwd', 'angular_flux_adj']:
        data[key] = np.transpose(data[key], (1, 2, 3, 4, 0))

    # Write pickles
    for key in sorted(data):
        #if key.startswith('angular'): continue
        write_pickle(data, key)

    return data

if __name__ == '__main__': main()
