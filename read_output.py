#!/usr/bin/python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
import itertools
import numpy as np
import h5py
from pyvisfile import silo

def main():
    data = read_all()
    print('')
    print('%-18s  %-7s  %-21s  %9s' % ('Key', 'Type', 'Shape', 'Size'))
    #sorted_keys = sorted(data)
    sorted_keys = sorted(sorted(data), key=lambda x: data[x].size, reverse=True)
    for key in sorted_keys:
        print('%-18s  %-7s  %-21s  %9s' %
              (key, data[key].dtype, data[key].shape, data[key].size))

def read_all():
    data = {}
    read_pickles(data)
    read_hdf5(data)
    read_silo(data)
    calculate(data)
    return data

def write_pickle(data, key):
    fname = 'pickles/%s.npy' % (key)
    print('Writing %s to file' % (fname))
    np.save(fname, data[key])

def read_pickles(data):
    print('Reading pickles from file')
    keys = ['source_indices', 'source_strengths', 'source_spectrum',
            'response_indices', 'response_strengths', 'response_spectrum',
            'mat_names', 'mix_table', 'material_map', 'sigma_t', 'sigma_s']
    for key in keys: data[key] = np.load('pickles/%s.npy' % (key))

def read_hdf5(data):
    keys_orig = data.keys()

    # Read HDF5 file containing forward flux
    fname = 'fwd_solution/denovo_forward.h5'
    print('Reading %s from file' % (fname))
    hf = h5py.File(fname,'r')

    # Load forward data
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
    keys_new = sorted([key for key in data.keys() if key not in keys_orig])
    for key in keys_new:
        #if key.startswith('angular'): continue
        write_pickle(data, key)

def read_silo(data):
    keys_orig = data.keys()

    # Read silo file
    fname = 'output/fields.silo'
    print('Reading %s from file' % (fname))
    db = silo.SiloFile(fname, create=False, mode=silo.DB_READ)
    toc = db.get_toc()

    # Dimensions
    nz  = len(data['mesh_z']) - 1  # Number of Z intervals
    ny  = len(data['mesh_y']) - 1  # Number of Y intervals
    nx  = len(data['mesh_x']) - 1  # Number of X intervals
    ngf = len(data['mesh_g'])      # Number of energy groups in flux

    # First energy group
    g0 = data['mesh_g'][0]

    # Initialize arrays
    data['flux_fwd' ] = np.zeros((nz, ny, nx, ngf))
    data['flux_adj' ] = np.zeros((nz, ny, nx, ngf))
    data['current_x'] = np.zeros((nz, ny, nx, ngf))
    data['current_y'] = np.zeros((nz, ny, nx, ngf))
    data['current_z'] = np.zeros((nz, ny, nx, ngf))

    # Read flux and current from file
    qvar_names = toc.qvar_names
    for name in qvar_names:
        key = None
        for i in data.keys():
            if name.startswith(i):
                key = i
                break
        if key is None: continue
        try: ig = int(name.split('_')[-1]) - g0
        except: continue
        vals = np.swapaxes(db.get_quadvar(name).vals[0], 0, 2)
        data[key][:, :, :, ig] = vals

    # Write pickles
    keys_new = sorted([key for key in data.keys() if key not in keys_orig])
    for key in keys_new: write_pickle(data, key)

def calculate(data):
    keys_orig = data.keys()

    # Dimensions
    n_mix = data['mix_table'       ].shape[0]  # Number of mixed materials
    nm    = data['mix_table'       ].shape[1]  # Number of pure materials
    nz    = data['angular_flux_fwd'].shape[0]  # Number of Z intervals
    ny    = data['angular_flux_fwd'].shape[1]  # Number of Y intervals
    nx    = data['angular_flux_fwd'].shape[2]  # Number of X intervals
    na    = data['angular_flux_fwd'].shape[3]  # Number of angles
    ngf   = data['angular_flux_fwd'].shape[4]  # Number of energy groups in flux
    ngx   = data['sigma_t'         ].shape[1]  # Number of energy groups in XS

    # First and last energy group
    g0 = data['mesh_g'][ 0]
    g1 = data['mesh_g'][-1]

    # Calculate source
    data['source'] = np.zeros((nz, ny, nx, ngf))
    for i, ind in enumerate(data['source_indices']):
        val = data['source_strengths'][i]
        iz = ind // (ny * nx)
        iy = (ind // nx) % ny
        ix = ind % nx
        data['source'][iz, iy, ix, :] = \
            val * data['source_spectrum'][g0:g1 + 1]

    # Calculate response
    data['response'] = np.zeros((nz, ny, nx, ngf))
    for i, ind in enumerate(data['response_indices']):
        val = data['response_strengths'][i]
        iz = ind // (ny * nx)
        iy = (ind // nx) % ny
        ix = ind % nx
        data['response'][iz, iy, ix, :] = \
            val * data['response_spectrum'][g0:g1 + 1]

    # Calculate contributon
    data['contributon'] = data['flux_fwd'] * data['flux_adj']

    # Calculate current magnitude
    data['current'] = np.sqrt((data['current_x']**2 + data['current_y']**2 +
                               data['current_z']**2))

    # Calculate total flux
    data['flux_fwd_int'   ] = np.sum(data['flux_fwd'   ], 3)
    data['flux_adj_int'   ] = np.sum(data['flux_adj'   ], 3)
    data['contributon_int'] = np.sum(data['contributon'], 3)

    # Calculate cross sections for mixed materials
    data['sigma_t_mixed'] = np.zeros((n_mix, ngx))
    data['sigma_s_mixed'] = np.zeros((n_mix, ngx, ngx))
    for i_mix, im in itertools.product(xrange(n_mix), xrange(nm)):
        frac = data['mix_table'][i_mix, im]
        if frac == 0.0: continue
        data['sigma_t_mixed'][i_mix, :] += data['sigma_t'][im, :] * frac
        data['sigma_s_mixed'][i_mix, :, :] += data['sigma_s'][im, :, :] * frac

    # Calculate perturbed cross sections
    data['sigma_t_pert'] = np.zeros((n_mix, nm, ngx))
    data['sigma_s_pert'] = np.zeros((n_mix, nm, ngx, ngx))
    for i_mix, im in itertools.product(xrange(n_mix), xrange(nm)):
        data['sigma_t_pert'][i_mix, im, :] = \
            data['sigma_t'][im, :] - data['sigma_t_mixed'][i_mix, :]
        data['sigma_s_pert'][i_mix, im, :, :] = \
            data['sigma_s'][im, :, :] - data['sigma_s_mixed'][i_mix, :, :]

    # Write pickles
    keys_new = sorted([key for key in data.keys() if key not in keys_orig])
    for key in keys_new: write_pickle(data, key)

if __name__ == '__main__': main()
