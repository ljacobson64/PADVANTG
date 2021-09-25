#!/usr/bin/python3

import os
import copy
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from cycler import cycler

from draw_geom import draw_geometry
from mcnp_colors import mcnp_colors as MC

small_plots = True
if small_plots:
    bbox_inches = 'tight'
    font_size = 16.0
else:
    bbox_inches = None
    font_size = 10.0

plt.rcParams['font.size'       ] = font_size
plt.rcParams['font.family'     ] = 'serif'
plt.rcParams['mathtext.default'] = 'regular'

mat_names_short = {
    'void'                            : 'Void'        ,
    'Air (Dry, Near Sea Level)'       : 'Air'         ,
    'Aluminum, Alloy 6061-O'          : 'Aluminum'    ,
    'Beryllium'                       : 'Beryllium'   ,
    'Beryllium Oxide'                 : 'BeO'         ,
    'Boron'                           : 'Boron'       ,
    'Carbon, Graphite (Reactor Grade)': 'Graphite'    ,
    'Concrete, Regular'               : 'Concrete'    ,
    'Copper'                          : 'Copper'      ,
    'Gadolinium'                      : 'Gadolinium'  ,
    'Iron'                            : 'Iron'        ,
    'Polyethylene, Borated'           : 'Borated HDPE',
    'Polyethylene, Non-borated'       : 'HDPE'        ,
    'Water, Heavy'                    : 'Heavy Water' ,
    'Water, Liquid'                   : 'Water'       ,
    'Zircaloy-4'                      : 'Zircaloy-4'  ,
    'Uranium-235'                     : 'Uranium-235' }

def main():
    data = read_hdf5()
    plot_all(data)

def read_hdf5():
    key_dict = {
        'fwd_solution/denovo-forward.inp.h5':
            {'group_bounds_n'    : 'group_bounds_n'           ,
             'group_bounds_p'    : 'group_bounds_p'           ,
             'mixtable'          : 'mixtable'                 ,
             'matids'            : 'matids'                   ,
             'source_spectra'    : 'volsrc/spectra'           ,
             'source_strength'   : 'volsrc/strength'          ,
             'source_ids'        : 'volsrc/ids'               },
        'adj_solution/denovo-adjoint.inp.h5':
            {'response_spectra'  : 'volsrc/spectra'           ,
             'response_strength' : 'volsrc/strength'          ,
             'response_ids'      : 'volsrc/ids'               },
        'fwd_solution/denovo-forward.out.h5':
            {'mesh_x'            : 'denovo/mesh_x'            ,
             'mesh_y'            : 'denovo/mesh_y'            ,
             'mesh_z'            : 'denovo/mesh_z'            ,
             'mesh_g'            : 'denovo/mesh_g'            ,
             'quadrature_angles' : 'denovo/quadrature_angles' ,
             'quadrature_weights': 'denovo/quadrature_weights'},
        'adj_solution/denovo-adjoint.out.h5': {},
        'custom_output/data1.h5'            :
            {'mat_names'         : 'mat_names'                ,
             'sigma_t'           : 'sigma_t'                  ,
             'sigma_s'           : 'sigma_s'                  },
        'custom_output/data3.h5'            :
            {'reverse_angle_map' : 'reverse_angle_map'        ,
             'source'            : 'source'                   ,
             'response'          : 'response'                 ,
             'scalar_flux_fwd'   : 'scalar_flux_fwd'          ,
             'scalar_flux_adj'   : 'scalar_flux_adj'          ,
             'scalar_flux_con'   : 'scalar_flux_con'          ,
             'sigma_t_mixed'     : 'sigma_t_mixed'            ,
             'sigma_s_mixed'     : 'sigma_s_mixed'            ,
             'sigma_t_pert'      : 'sigma_t_pert'             ,
             'sigma_s_pert'      : 'sigma_s_pert'             ,
             'current_fwd'       : 'current_fwd'              ,
             'current_adj'       : 'current_adj'              ,
             'current_con'       : 'current_con'              ,
             'dR'                : 'dR'                       }
    }

    data = {}
    for fname in key_dict.keys():
        hf = h5py.File(fname, 'r')
        for key, key_h5 in key_dict[fname].items(): data[key] = hf[key_h5][()]
    return data

def plot_all(data):
    nm     = data['dR'].shape[0]  # Number of pure materials
    nz     = data['dR'].shape[1]  # Number of Z intervals
    g0     = data['mesh_g'][0]    # First energy group
    x_vals = data['mesh_x']       # Horizontal direction
    y_vals = data['mesh_y']       # Vertical direction
    hz     = nz // 2              # Use the middle Z slice

    # Plot spectra
    plot_spectra(data['source_spectra'], data['response_spectra'])

    # All future drawn lines should be black
    plt.rcParams['axes.prop_cycle'] = cycler(color=['k'])

    # Plot material map
    plot_data = data['matids'][hz, :, :]
    plot_mask = np.where(plot_data != 14)
    plot_map(plot_data,
             x_vals, y_vals,
             'Material Map',
             'material_map.png',
             'mats',
             mat_names=data['mat_names'])

    # Plot total source
    plot_map(np.sum(data['source'][:, hz, :, :], 0),
             x_vals, y_vals,
             'Source (Total)',
             'source_total.png',
             'lin')

    # Plot total response
    plot_map(np.sum(data['response'][:, hz, :, :], 0),
             x_vals, y_vals,
             'Response (Total)',
             'response_total.png',
             'lin')

    # Plot total forward flux
    plot_map(np.sum(data['scalar_flux_fwd'][:, hz, :, :], 0),
             x_vals, y_vals,
             'Forward Flux (Total)',
             'scalar_flux_fwd_total.png',
             'log')

    # Plot total adjoint flux
    plot_map(np.sum(data['scalar_flux_adj'][:, hz, :, :], 0),
             x_vals, y_vals,
             'Adjoint Flux (Total)',
             'scalar_flux_adj_total.png',
             'log')

    # Plot total contributon flux
    plot_map(np.sum(data['scalar_flux_con'][:, hz, :, :], 0),
             x_vals, y_vals,
             'Contributon Flux (Total)',
             'scalar_flux_con_total.png',
             'log')

    # Plot total forward current
    plot_quiver(np.sum(data['current_fwd'][:, hz, :, :, :2], 0),
                x_vals, y_vals,
                'Forward Current (Total)',
                'current_fwd_total.png')

    # Plot total adjoint current
    plot_quiver(np.sum(data['current_adj'][:, hz, :, :, :2], 0),
                x_vals, y_vals,
                'Adjoint Current (Total)',
                'current_adj_total.png')

    # Plot total contributon current
    plot_quiver(np.sum(data['current_con'][:, hz, :, :, :2], 0),
                x_vals, y_vals,
                'Contributon Current (Total)',
                'current_con_total.png')

    # Plots for the fast group (2) and thermal group (25)
    for igx in [2, 25]:
        igf = igx - g0

        # Forward flux
        plot_map(data['scalar_flux_fwd'][igf, hz, :, :],
                 x_vals, y_vals,
                 'Forward Flux (Group %u)'   % (igx),
                 'scalar_flux_fwd_g%02u.png' % (igx),
                 'log')

        # Adjoint flux
        plot_map(data['scalar_flux_adj'][igf, hz, :, :],
                 x_vals, y_vals,
                 'Adjoint Flux (Group %u)'   % (igx),
                 'scalar_flux_adj_g%02u.png' % (igx),
                 'log')

        # Contributon flux
        plot_map(data['scalar_flux_con'][igf, hz, :, :],
                 x_vals, y_vals,
                 'Contributon Flux (Group %u)' % (igx),
                 'scalar_flux_con_g%02u.png'   % (igx),
                 'log')

        # Forward current
        plot_quiver(data['current_fwd'][igf, hz, :, :, :2],
                    x_vals, y_vals,
                    'Forward Current (Group %u)' % (igx),
                    'current_fwd_g%02u.png'      % (igx))

        # Adjoint current
        plot_quiver(data['current_adj'][igf, hz, :, :, :2],
                    x_vals, y_vals,
                    'Adjoint Current (Group %u)' % (igx),
                    'current_adj_g%02u.png'      % (igx))

        # Contributon current
        plot_quiver(data['current_con'][igf, hz, :, :, :2],
                    x_vals, y_vals,
                    'Contributon Current (Group %u)' % (igx),
                    'current_con_g%02u.png'          % (igx))

    # Plot dR for all materials
    for im in range(nm):
        mat_name_short = mat_names_short[data['mat_names'][im].decode()]
        plot_map(data['dR'][im, hz, :, :],
                 x_vals, y_vals,
                 r'$\delta R$ for %s' % (mat_name_short),
                 'dR_%02u.png' % (im),
                 'logplusminus',
                 plot_mask=plot_mask)

def plot_map(plot_data, x_vals, y_vals, title, fname, fmt, plot_mask=None,
             mat_names=None):
    # Make a copy of the plot data to avoid overwriting it
    plot_data_use = np.array(plot_data)

    # Apply mask if appropriate
    if plot_mask is not None: plot_data_use[plot_mask] = 0.0

    # Color scale
    if fmt == 'logplusminus':
        logvmax = 8
        span = 8
        loglinthresh = logvmax - span
        vmax = 10**logvmax
        vmin = -vmax
        linthresh = 10**loglinthresh
        norm = colors.SymLogNorm(linthresh=linthresh, linscale=1.0,
                                 vmin=vmin, vmax=vmax, base=10.0)
        cmap='RdBu_r'
        ticks_neg = [-10**x for x in
                     np.linspace(loglinthresh, logvmax, span + 1)]
        ticks_pos = [10**x for x in
                     np.linspace(loglinthresh, logvmax, span + 1)]
        ticks = list(reversed(ticks_neg)) + [0.0] + ticks_pos
        tick_labels_neg = [r'$-10^{%u}$' % (x) for x in
                           np.arange(loglinthresh, logvmax + 1, 1)]
        tick_labels_pos = [r'$10^{%u}$'  % (x) for x in
                           np.arange(loglinthresh, logvmax + 1, 1)]
        tick_labels = (list(reversed(tick_labels_neg)) + [0] + tick_labels_pos)
        #print('Min/max values: %9.2e, %9.2e' %
        #      (np.min(plot_data_use), np.max(plot_data_use)))
    elif fmt == 'log':
        result_min_nonzero = np.min(plot_data_use[plot_data_use != 0])
        result_max_nonzero = np.max(plot_data_use[plot_data_use != 0])
        logvmax = np.ceil(np.log10(result_max_nonzero))
        logvmin = max(np.floor(np.log10(result_min_nonzero)), logvmax - 10)
        vmin = 10**logvmin
        vmax = 10**logvmax
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        cmap = copy.copy(cm.get_cmap('Spectral_r'))
        cmap.set_bad(cmap(0.0))
        ticks = [10**x for x in np.arange(logvmin, logvmax + 1, 1)]
        tick_labels = [r'$10^{%u}$' % (x) for x in
                       np.arange(logvmin, logvmax + 1, 1)]
    elif fmt == 'lin':
        vmin = 0.0
        vmax = np.max(plot_data_use) * 1.1
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = 'Reds'
    elif fmt == 'mats':
        with open('com_1', 'r') as reader: lines = reader.readlines()
        unique_vals = np.unique(plot_data_use)
        for i, unique_val in enumerate(unique_vals):
            plot_data_use[np.where(plot_data_use == unique_val)] = i
        color_strs = (['white'] +
                      [x.split()[2] for x in lines if x.startswith('shade')])
        color_strs = [x for i, x in enumerate(color_strs) if i in unique_vals]
        CL = [MC[x] for x in color_strs]
        vmin = 0
        vmax = len(CL)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(CL)
        ticks = np.linspace(0.5, vmax - 0.5, len(CL))
        tick_labels = [mat_names_short[x.decode()]
                       for x in mat_names[unique_vals]]

    fig, ax = plt.subplots(1, 1)
    if 'mats' in fmt:
        fig.set_size_inches(10, 6)
        plt.subplots_adjust(left=-0.025, right=0.75)
    else:
        fig.set_size_inches(8, 6)
        plt.subplots_adjust(left=0.125, right=0.9)

    # Plot data
    im = ax.pcolormesh(x_vals, y_vals, plot_data_use, norm=norm, cmap=cmap)

    # Formatting
    ax.set_aspect('equal', 'box')
    ax.set_xlim([x_vals[0], x_vals[-1]])
    ax.set_ylim([y_vals[0], y_vals[-1]])
    ax.set_xticks(np.linspace(x_vals[0], x_vals[-1], 5))
    ax.set_yticks(np.linspace(y_vals[0], y_vals[-1], 5))
    ax.set_xlabel('x position [cm]')
    ax.set_ylabel('y position [cm]')
    ax.set_title(title)

    # Colorbar
    cbar = plt.colorbar(im)
    if 'log' in fmt:
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels(tick_labels)
    elif 'mats' in fmt:
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels(tick_labels)
        cbar.ax.tick_params(length=0)
        cbar.ax.invert_yaxis()

    # Draw geometry
    draw_geometry()

    # Save figure
    if not os.path.exists('images'): os.makedirs('images')
    print('Writing images/%s to file' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)
    plt.close()

def plot_quiver(plot_data, x_vals, y_vals, title, fname):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)
    plt.subplots_adjust(left=0.125, right=0.9)

    # Plot data
    x_centers = 0.5 * (x_vals[:-1] + x_vals[1:])
    y_centers = 0.5 * (y_vals[:-1] + y_vals[1:])
    im = ax.quiver(x_centers, y_centers, plot_data[:, :, 0], plot_data[:, :, 1])

    # Formatting
    ax.set_aspect('equal', 'box')
    ax.set_xlim([x_vals[0], x_vals[-1]])
    ax.set_ylim([y_vals[0], y_vals[-1]])
    ax.set_xticks(np.linspace(x_vals[0], x_vals[-1], 5))
    ax.set_yticks(np.linspace(y_vals[0], y_vals[-1], 5))
    ax.set_xlabel('x position [cm]')
    ax.set_ylabel('y position [cm]')
    ax.set_title(title)

    # Draw geometry
    draw_geometry()

    # Save figure
    if not os.path.exists('images'): os.makedirs('images')
    print('Writing images/%s to file' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)
    plt.close()

def plot_spectra(source, response):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    # Plot spectra
    ng = source.shape[1]
    ebins = range(ng)
    plt.bar(ebins, source  [0, :], label='Source'  )
    plt.bar(ebins, response[0, :], label='Response')

    # Formatting
    ax.set_xlim([-0.6, ng - 0.4])
    ax.set_xlabel('Energy group')
    ax.set_ylabel('Magnitude')
    ax.set_title('Source and Response Energy Spectra')
    ax.grid()
    ax.legend()

    # Save figure
    if not os.path.exists('images'): os.makedirs('images')
    fname = 'spectra_lin.png'
    print('Writing images/%s to file' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)

    # Log scale
    ax.set_yscale('log')

    # Save figure
    fname = 'spectra_log.png'
    print('Writing images/%s to file' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)

    plt.close()

if __name__ == '__main__': main()
