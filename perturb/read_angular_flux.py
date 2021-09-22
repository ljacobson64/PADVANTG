#!/usr/bin/python3

import os
import time
import numpy as np
import h5py

def add_array(hf, name, data):
    hf.create_dataset(name, data=data, compression='gzip', compression_opts=0)

def main():
    # Read HDF5 file containing forward flux
    time_start = time.time()
    fname = 'fwd_solution/denovo_forward.h5'
    print('Reading %s from file...' % (fname), end='', flush=True)
    hf = h5py.File(fname, 'r')

    # Load forward data
    data = {}
    for key in hf:
        if key in ['angles', 'quadrature_weights'] or key.startswith('mesh'):
            new_key = str(key)
        elif key in ['angular_flux']: new_key = 'angular_flux_fwd'
        else: continue
        data[new_key] = hf.get(key)[()]
    hf.close()
    time_end = time.time()
    elapsed = time_end - time_start
    print(' took %7.3f seconds' % (elapsed))

    # Read HDF5 file containing adjoint flux
    time_start = time.time()
    fname = 'adj_solution/denovo_adjoint.h5'
    print('Reading %s from file...' % (fname), end='', flush=True)
    hf = h5py.File(fname, 'r')

    # Load adjoint data
    data['angular_flux_adj'] = hf.get('angular_flux')[()]
    hf.close()
    time_end = time.time()
    elapsed = time_end - time_start
    print(' took %7.3f seconds' % (elapsed))

    # Move energy axis to end
    for key in ['angular_flux_fwd', 'angular_flux_adj']:
        data[key] = np.transpose(data[key], (1, 2, 3, 4, 0))

    # Write arrays to HDF5
    time_start = time.time()
    if not os.path.exists('custom_output'): os.makedirs('custom_output')
    fname = 'custom_output/data2.h5'
    print('Writing data to %s...          ' % (fname), end='', flush=True)
    hf = h5py.File(fname, 'w')
    for key in sorted(data):
        #if key.startswith('angular'): continue
        add_array(hf, key, data[key])
    hf.close()
    time_end = time.time()
    elapsed = time_end - time_start
    print(' took %7.3f seconds' % (elapsed))

if __name__ == '__main__': main()
