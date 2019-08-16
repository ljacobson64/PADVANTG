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
    sorted_keys = sorted(data)
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

    # Change axes of angular flux
    for key in ['angular_flux_fwd', 'angular_flux_adj']:
        data[key] = np.transpose(data[key], (3, 2, 1, 4, 0))

    # Write pickles
    keys_final = data.keys()
    keys_new = [key for key in keys_final if key not in keys_orig]
    for key in keys_new:
        #if key.startswith('angular'): continue
        np.save('pickles/%s' % (key), data[key])

def read_silo(data):
    keys_orig = data.keys()

    # Read silo file
    fname = 'output/fields.silo'
    print('Reading %s from file' % (fname))
    db = silo.SiloFile(fname, create=False, mode=silo.DB_READ)
    toc = db.get_toc()

    # Mesh
    ng = len(data['mesh_g'])
    nx = len(data['mesh_x']) - 1
    ny = len(data['mesh_y']) - 1
    nz = len(data['mesh_z']) - 1
    g0 = data['mesh_g'][0]

    # Initialize arrays
    data['flux_fwd' ] = np.zeros((ng, nx, ny, nz))
    data['flux_adj' ] = np.zeros((ng, nx, ny, nz))
    data['current_x'] = np.zeros((ng, nx, ny, nz))
    data['current_y'] = np.zeros((ng, nx, ny, nz))
    data['current_z'] = np.zeros((ng, nx, ny, nz))

    # Read flux and current from file
    qvar_names = toc.qvar_names
    for name in qvar_names:
        key = None
        for i in data.keys():
            if name.startswith(i):
                key = i
                break
        if key is None: continue
        try: group = int(name.split('_')[-1]) - g0
        except: group = None
        vals = db.get_quadvar(name).vals[0]
        if group is None: data[key] = vals
        else: data[key][group, :, :, :] = vals

    # Write pickles
    keys_final = data.keys()
    keys_new = [key for key in keys_final if key not in keys_orig]
    for key in keys_new: np.save('pickles/%s' % (key), data[key])

def calculate(data):
    keys_orig = data.keys()

    ng_tot    = data['sigma_t'  ].shape[1]
    ng_flux   = data['flux_fwd' ].shape[0]
    nx        = data['flux_fwd' ].shape[1]
    ny        = data['flux_fwd' ].shape[2]
    nz        = data['flux_fwd' ].shape[3]
    num_mixed = data['mix_table'].shape[0]
    num_mats  = data['mix_table'].shape[1]

    g0 = data['mesh_g'][0]
    g1 = data['mesh_g'][-1]

    # Calculate source
    data['source'] = np.zeros((ng_flux, nx, ny, nz))
    for i, ind in enumerate(data['source_indices']):
        val = data['source_strengths'][i]
        ix = ind % nx
        iy = (ind // nx) % ny
        iz = ind // (ny * nx)
        data['source'][:, ix, iy, iz] = \
            val * data['source_spectrum'][g0:g1 + 1]

    # Calculate response
    data['response'] = np.zeros((ng_flux, nx, ny, nz))
    for i, ind in enumerate(data['response_indices']):
        val = data['response_strengths'][i]
        ix = ind % nx
        iy = (ind // nx) % ny
        iz = ind // (ny * nx)
        data['response'][:, ix, iy, iz] = \
            val * data['response_spectrum'][g0:g1 + 1]

    # Calculate contributon
    data['contributon'] = data['flux_fwd'] * data['flux_adj']

    # Calculate current magnitude
    data['current'] = np.sqrt((data['current_x']**2 + data['current_y']**2 +
                               data['current_z']**2))

    # Calculate total flux
    data['flux_fwd_int'   ] = np.sum(data['flux_fwd'   ], 0)
    data['flux_adj_int'   ] = np.sum(data['flux_adj'   ], 0)
    data['contributon_int'] = np.sum(data['contributon'], 0)

    # Calculate cross sections for mixed materials
    data['sigma_t_mixed'] = np.zeros((num_mixed, ng_tot))
    data['sigma_s_mixed'] = np.zeros((num_mixed, ng_tot, ng_tot))
    for i_mix, i_mat in itertools.product(xrange(num_mixed), xrange(num_mats)):
        frac = data['mix_table'][i_mix, i_mat]
        if frac == 0.0: continue
        data['sigma_t_mixed'][i_mix, :] += \
            data['sigma_t'][i_mat, :] * frac
        data['sigma_s_mixed'][i_mix, :, :] += \
            data['sigma_s'][i_mat, :, :] * frac

    # Calculate perturbed cross sections
    data['sigma_t_pert'] = np.zeros((num_mixed, num_mats, ng_tot))
    data['sigma_s_pert'] = np.zeros((num_mixed, num_mats, ng_tot, ng_tot))
    for i_mix, i_mat in itertools.product(xrange(num_mixed), xrange(num_mats)):
        data['sigma_t_pert'][i_mix, i_mat, :] = \
            data['sigma_t'][i_mat, :] - data['sigma_t_mixed'][i_mix, :]
        data['sigma_s_pert'][i_mix, i_mat, :, :] = \
            data['sigma_s'][i_mat, :, :] - data['sigma_s_mixed'][i_mix, :, :]

    # Write pickles
    keys_final = data.keys()
    keys_new = [key for key in keys_final if key not in keys_orig]
    for key in keys_new: np.save('pickles/%s' % (key), data[key])

if __name__ == '__main__': main()
