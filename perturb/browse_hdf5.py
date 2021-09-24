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
    fnames = ['fwd_solution/denovo-forward.inp.h5',
              'adj_solution/denovo-adjoint.inp.h5',
              'fwd_solution/denovo-forward.out.h5',
              'adj_solution/denovo-adjoint.out.h5']

    hriz_line = ('+%s+%s+%s+%s+%s+%s+' %
                 ('-' * 23, '-' * 27, '-' * 23, '-' * 9, '-' * 13, '-' * 27))
    print(hriz_line)
    print('| %-21s | %-25s | %-21s | %-7s | %11s | %-25s |' %
          ('Filename', 'Key', 'Shape', 'Type', 'Num_Nonzero', 'Label'))
    print(hriz_line)

    for fname in fnames:
        hf = h5py.File(fname, 'r')
        data = load_hdf5_data(hf)

        if 'mixtable' in data.keys():
            data['mixtable'] = np.array([[float(x) for x in data['mixtable'][i]]
                                        for i in range(len(data['mixtable']))])

        labels = OrderedDict()
        if fname == 'fwd_solution/denovo-forward.inp.h5':
            labels['group_bounds_n'           ] = 'Energy groups (neutron)'
            labels['group_bounds_p'           ] = 'Energy groups (photon)'
            labels['mixtable'                 ] = 'Mix table'
            labels['matids'                   ] = 'Material map'
            labels['volsrc/spectra'           ] = 'Source spectra'
            labels['volsrc/strength'          ] = 'Source strength'
            labels['volsrc/ids'               ] = 'Source IDs'
        elif fname == 'adj_solution/denovo-adjoint.inp.h5':
            labels['volsrc/spectra'           ] = 'Response spectra'
            labels['volsrc/strength'          ] = 'Response strength'
            labels['volsrc/ids'               ] = 'Response IDs'
        elif fname == 'fwd_solution/denovo-forward.out.h5':
            labels['denovo/mesh_x'            ] = 'Mesh points (x)'
            labels['denovo/mesh_y'            ] = 'Mesh points (y)'
            labels['denovo/mesh_z'            ] = 'Mesh points (z)'
            labels['denovo/mesh_g'            ] = 'Mesh points (energy)'
            labels['denovo/quadrature_angles' ] = 'Quadrature angles'
            labels['denovo/quadrature_weights'] = 'Quadrature weights'
            labels['denovo/angular_flux'      ] = 'Angular flux (forward)'
        elif fname == 'adj_solution/denovo-adjoint.out.h5':
            labels['denovo/angular_flux'      ] = 'Angular flux (adjoint)'

        for key, label in labels.items():
            key_list = key.split('/')
            val = data[key_list[0]]
            for key2 in key_list[1:]: val = val[key2]
            num_nonzero = np.count_nonzero(val)
            print('| %-21s | %-25s | %-21s | %-7s | %11s | %-25s |' %
                  (fname.split('/')[-1], key, val.shape, val.dtype,
                   num_nonzero, label))
        print(hriz_line)

if __name__ == '__main__': main()
