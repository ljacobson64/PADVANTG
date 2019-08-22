#!/usr/bin/python

import itertools
import numpy as np

# Flag to calculate contributon with angular flux
use_angular = True

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
    calculate(data)
    return data

def write_pickle(data, key):
    fname = 'pickles/%s.npy' % (key)
    print('Writing %s to file' % (fname))
    np.save(fname, data[key])

def read_pickles(data):
    keys = ['source_indices'  , 'source_spectrum'  , 'source_strengths'  ,
            'response_indices', 'response_spectrum', 'response_strengths',
            'mesh_x', 'mesh_y', 'mesh_z', 'mesh_g', 'angles',
            'mat_names', 'material_map', 'mix_table', 'sigma_t', 'sigma_s',
            'quadrature_weights', 'angular_flux_fwd', 'angular_flux_adj']
    for key in keys:
        fname = 'pickles/%s.npy' % (key)
        print('Reading %s from file' % (fname))
        data[key] = np.load(fname)

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

    # Calculate reverse angle map
    data['reverse_angle_map'] = np.zeros(na, dtype=np.int64)
    for i in xrange(na):
        ix = data['angles'][i, 0]
        iy = data['angles'][i, 1]
        iz = data['angles'][i, 2]
        for j in xrange(na):
            jx = data['angles'][j, 0]
            jy = data['angles'][j, 1]
            jz = data['angles'][j, 2]
            if ix == -jx and iy == -jy and iz == -jz:
                data['reverse_angle_map'][i] = j
                break

    # Calculate source
    print('Calculating source')
    data['source'] = np.zeros((nz, ny, nx, ngf))
    for i, ind in enumerate(data['source_indices']):
        val = data['source_strengths'][i]
        iz = ind // (ny * nx)
        iy = (ind // nx) % ny
        ix = ind % nx
        data['source'][iz, iy, ix, :] = \
            val * data['source_spectrum'][g0:g1 + 1]

    # Calculate response
    print('Calculating response')
    data['response'] = np.zeros((nz, ny, nx, ngf))
    for i, ind in enumerate(data['response_indices']):
        val = data['response_strengths'][i]
        iz = ind // (ny * nx)
        iy = (ind // nx) % ny
        ix = ind % nx
        data['response'][iz, iy, ix, :] = \
            val * data['response_spectrum'][g0:g1 + 1]

    # Calculate scalar flux
    print('Calculating scalar flux')
    qw = data['quadrature_weights'][None, None, None, :, None] / (4.0 * np.pi)
    data['flux_fwd'] = np.sum(data['angular_flux_fwd'] * qw, 3)
    data['flux_adj'] = np.sum(data['angular_flux_adj'] * qw, 3)

    # Calculate contributon
    print('Calculating contributon')
    if use_angular:
        data['contributon'] = np.sum(data['angular_flux_fwd'] *
            data['angular_flux_adj'][:, :, :, data['reverse_angle_map'], :] *
            qw, 3)
    else: data['contributon'] = data['flux_fwd'] * data['flux_adj']

    # Calculate energy-integrated flux
    print('Calculating energy-integrated flux')
    data['flux_fwd_int'   ] = np.sum(data['flux_fwd'   ], 3)
    data['flux_adj_int'   ] = np.sum(data['flux_adj'   ], 3)
    data['contributon_int'] = np.sum(data['contributon'], 3)

    # Calculate cross sections for mixed materials
    print('Calculating cross sections for mixed materials')
    data['sigma_t_mixed'] = np.zeros((n_mix, ngx))
    data['sigma_s_mixed'] = np.zeros((n_mix, ngx, ngx))
    for i_mix, im in itertools.product(xrange(n_mix), xrange(nm)):
        frac = data['mix_table'][i_mix, im]
        if frac == 0.0: continue
        data['sigma_t_mixed'][i_mix, :] += data['sigma_t'][im, :] * frac
        data['sigma_s_mixed'][i_mix, :, :] += data['sigma_s'][im, :, :] * frac

    # Calculate perturbed cross sections
    print('Calculating perturbed cross sections')
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
