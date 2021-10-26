#!/usr/bin/python3

import os
import sys
import copy
from collections import OrderedDict
import argparse
import multiprocessing as mp
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from cycler import cycler


from mcnp_colors import get_mcnp_color

sys.path.insert(0, os.getcwd())
from draw_geom import draw_geometry

small_plots = False
if small_plots:
    bbox_inches = 'tight'
    font_size = 16.0
else:
    bbox_inches = None
    font_size = 10.0

plt.rcParams['font.size'       ] = font_size
plt.rcParams['font.family'     ] = 'serif'
plt.rcParams['mathtext.default'] = 'regular'

def main():
    data = read_hdf5()
    plot_args_list = get_all_plot_args(data)

    jobs = get_args().jobs
    if jobs > 1: pool = mp.Pool(processes=jobs)
    for plot_args in plot_args_list:
        if jobs > 1: pool.apply_async(make_plot, args=(plot_args,))
        else: make_plot(plot_args)
    if jobs > 1:
        pool.close()
        pool.join()

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help='number of CPUs')
    args = parser.parse_args()
    return args

def add_to_key_dict(key_dict, tokens):
    fname       = tokens[0]
    key_hdf5    = tokens[1]
    key_import  = tokens[2]
    if fname in key_dict.keys(): pass
    else: key_dict[fname] = []
    key_dict[fname].append((key_hdf5, key_import))

def read_hdf5():
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

    # Load datasets
    data = {}
    for fname in key_dict.keys():
        print('Reading %s' % (fname))
        hf = h5py.File(fname, 'r')
        for key_hdf5, key_import in key_dict[fname]:
            if key_import.startswith('angular'): continue
            if '[' in key_import and ']' in key_import:
                ki1 =     key_import.split('[')[0].strip()
                ki2 = int(key_import.split('[')[1].split(']')[0].strip())
                if ki1 in data.keys(): data[ki1][ki2] = hf[key_hdf5][()]
                else: data[ki1] = {ki2: hf[key_hdf5][()]}
            else: data[key_import] = hf[key_hdf5][()]
    return data

def make_plot(plot_args):
    func = plot_args[0]
    args = plot_args[1:]
    func(*args)

def get_all_plot_args(data):
    nm     = data['mat_names'].shape[0]   # Number of pure materials
    nz     = data['mesh_z'].shape[0] - 1  # Number of Z intervals
    ng     = data['sigma_t'].shape[1]     # Number of energy bins
    x_vals = data['mesh_x']               # Horizontal direction
    y_vals = data['mesh_y']               # Vertical direction
    hz     = nz // 2                      # Use the middle Z slice

    # Get tally IDs
    tally_list = sorted(data['source_spectra_adj'].keys())

    # Get tally ratio IDs
    if os.path.exists('tally_ratio_list.txt'):
        with open('tally_ratio_list.txt', 'r') as r:
            tallies_ratio_list = [[int(x) for x in x.split()]
                                  for x in r.readlines()]
    else: tallies_ratio_list = []

    # List of all plot arguments
    plots = []

    # Plot spectra
    plots.append((plot_spectra,
                  data['source_spectra_fwd'],
                  data['source_spectra_adj']))

    # Plot material map
    plots.append((plot_map,
                  data['matids'][hz, :, :],
                  x_vals, y_vals,
                  'Material Map',
                  'material_map.png',
                  'mats',
                  None, None,
                  data['mat_names']))

    # Plot total forward source
    plots.append((plot_map,
                  np.sum(data['source_fwd'][:, hz, :, :], 0),
                  x_vals, y_vals,
                  'Forward Source (Total)',
                  'source_fwd/total.png',
                  'lin'))

    # Plot total adjoint response
    for ias, tally in enumerate(tally_list):
        plots.append((plot_map,
                      np.sum(data['source_adj'][ias, :, hz, :, :], 0),
                      x_vals, y_vals,
                      'Adjoint Source (Tally %u, Total)' % (tally),
                      'source_adj/t%03u_total.png'       % (tally),
                      'lin'))

    # Plot total forward flux
        plots.append((plot_map,
                      np.sum(data['scalar_flux_fwd'][:, hz, :, :], 0),
                      x_vals, y_vals,
                      'Forward Flux (Total)',
                      'scalar_flux_fwd/total.png',
                      'log'))

    # Plot total adjoint flux
    for ias, tally in enumerate(tally_list):
        plots.append((plot_map,
                      np.sum(data['scalar_flux_adj'][ias, :, hz, :, :], 0),
                      x_vals, y_vals,
                      'Adjoint Flux (Tally %u, Total)'  % (tally),
                      'scalar_flux_adj/t%03u/total.png' % (tally),
                      'log'))

    # Plot total contributon flux
    for ias, tally in enumerate(tally_list):
        plots.append((plot_map,
                      np.sum(data['scalar_flux_con'][ias, :, hz, :, :], 0),
                      x_vals, y_vals,
                      'Contributon Flux (Tally %u, Total)' % (tally),
                      'scalar_flux_con/t%03u/total.png'    % (tally),
                      'log'))

    # Plot total forward current
    plots.append((plot_quiver,
                  np.sum(data['current_fwd'][:, hz, :, :, :2], 0),
                  x_vals, y_vals,
                  'Forward Current (Total)',
                  'current_fwd/total.png'))

    # Plot total adjoint current
    for ias, tally in enumerate(tally_list):
        plots.append((plot_quiver,
                      np.sum(data['current_adj'][ias, :, hz, :, :, :2], 0),
                      x_vals, y_vals,
                      'Adjoint Current (Tally %u, Total)' % (tally),
                      'current_adj/t%03u/total.png'       % (tally)))

    # Plot total contributon current
    for ias, tally in enumerate(tally_list):
        plots.append((plot_quiver,
                      np.sum(data['current_con'][ias, :, hz, :, :, :2], 0),
                      x_vals, y_vals,
                      'Contributon Current (Tally %u, Total)' % (tally),
                      'current_con/t%03u/total.png'           % (tally)))

    # Plots for all energy groups
    for igx in range(ng):
        # Forward flux
        plots.append((plot_map,
                      data['scalar_flux_fwd'][igx, hz, :, :],
                      x_vals, y_vals,
                      'Forward Flux (Group %u)'   % (igx),
                      'scalar_flux_fwd/g%03u.png' % (igx),
                      'log',
                      np.max(data['scalar_flux_fwd'][:, hz, :, :])))

        # Adjoint flux
        for ias, tally in enumerate(tally_list):
            plots.append((plot_map,
                          data['scalar_flux_adj'][ias, igx, hz, :, :],
                          x_vals, y_vals,
                          'Adjoint Flux (Tally %u, Group %u)' % (tally, igx),
                          'scalar_flux_adj/t%03u/g%03u.png'  % (tally, igx),
                          'log',
                          np.max(data['scalar_flux_adj'][ias, :, hz, :, :])))

        # Contributon flux
        for ias, tally in enumerate(tally_list):
            plots.append((plot_map,
                          data['scalar_flux_con'][ias, igx, hz, :, :],
                          x_vals, y_vals,
                          'Contributon Flux (Tally %u, Group %u)' % (tally, igx),
                          'scalar_flux_con/t%03u/g%03u.png'       % (tally, igx),
                          'log',
                          np.max(data['scalar_flux_con'][ias, :, hz, :, :])))

        # Forward current
        plots.append((plot_quiver,
                      data['current_fwd'][igx, hz, :, :, :2],
                      x_vals, y_vals,
                      'Forward Current (Group %u)' % (igx),
                      'current_fwd/g%03u.png'      % (igx)))

        # Adjoint current
        for ias, tally in enumerate(tally_list):
            plots.append((plot_quiver,
                          data['current_adj'][ias, igx, hz, :, :, :2],
                          x_vals, y_vals,
                          'Adjoint Current (Tally %u, Group %u)' % (tally, igx),
                          'current_adj/t%03u/g%03u.png'          % (tally, igx)))

        # Contributon current
        for ias, tally in enumerate(tally_list):
            plots.append((plot_quiver,
                          data['current_con'][ias, igx, hz, :, :, :2],
                          x_vals, y_vals,
                          'Contributon Current (Tally %u, Group %u)' % (tally, igx),
                          'current_con/t%03u/g%03u.png'              % (tally, igx)))

    # Plot dR for all materials
    for ias, tally in enumerate(tally_list):
        for im in range(nm):
            mat_name = get_mat_name_short(data['mat_names'][im].decode())
            plots.append((plot_map,
                          data['dR'][ias, im, hz, :, :],
                          x_vals, y_vals,
                          r'$\delta R$ (Tally %u, %s)' % (tally, mat_name),
                          'dR/t%03u/m%03u.png'         % (tally, im      ),
                          'logplusminus',
                          np.max(np.abs(data['dR'][ias, :, hz, :, :]))))

    # Plot dR ratios for all materials
    for tally1, tally2 in tallies_ratio_list:
        ias1 = tally_list.index(tally1)
        ias2 = tally_list.index(tally2)
        ratio_data = (data['dR'][ias1, :, hz, :, :] /
                      data['dR'][ias2, :, hz, :, :])
        ratio_data[np.isnan(ratio_data)] = 0.0
        vmax = np.max(np.abs(ratio_data))
        for im in range(nm):
            mat_name = get_mat_name_short(data['mat_names'][im].decode())
            plot_data = ratio_data[im, :, :]
            plots.append((plot_map,
                          plot_data,
                          x_vals, y_vals,
                          r'$\delta R$ (Tally Ratio %u/%u, %s)' % (tally1, tally2, mat_name),
                          'dR_ratios/t%03u_t%03u/m%03u.png'     % (tally1, tally2, im      ),
                          'logplusminus',
                          vmax))

    return plots

def plot_map(plot_data, x_vals, y_vals, title, fname, fmt,
             vmax=None, plot_mask=None, mat_names=None):
    # Make a copy of the plot data to avoid overwriting it
    plot_data_use = np.array(plot_data)

    # Apply mask if appropriate
    if plot_mask is not None: plot_data_use[plot_mask] = 0.0

    # Color scale
    if fmt == 'logplusminus':
        if vmax is None: logvmax = 8
        else:            logvmax = int(np.ceil(np.log10(vmax)))
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
        if vmax is None:
            result_min_nonzero = np.min(plot_data_use[plot_data_use != 0])
            result_max_nonzero = np.max(plot_data_use[plot_data_use != 0])
            logvmax = np.ceil(np.log10(result_max_nonzero))
        else:
            logvmax = np.ceil(np.log10(vmax))
        logvmin = logvmax - 10
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
        if vmax is None: vmax = np.max(plot_data_use) * 1.1
        else:            vmax = vmax * 1.1
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = 'Reds'
    elif fmt == 'mats':
        with open('com_z', 'r') as reader: lines = reader.readlines()
        unique_vals = np.unique(plot_data_use)
        for i, unique_val in enumerate(unique_vals):
            plot_data_use[np.where(plot_data_use == unique_val)] = i
        color_strs = (['white'] +
                      [x.split()[2] for x in lines if x.startswith('shade')])
        color_strs = [x for i, x in enumerate(color_strs) if i in unique_vals]
        tick_labels = [get_mat_name_short(x.decode()) for x in
                       mat_names[unique_vals[:len(color_strs)]]]
        if np.sum(plot_data_use > len(color_strs) - 1) > 0:
            color_strs.append('white')
            tick_labels.append('(Mix)')
        CL = [get_mcnp_color(x) for x in color_strs]
        vmin = 0
        vmax = len(CL)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(CL)
        ticks = np.linspace(0.5, vmax - 0.5, len(CL))

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
    if not os.path.exists('images'): os.system('mkdir -p images')
    if '/' in fname:
        os.system('mkdir -p images/%s' % ('/'.join(fname.split('/')[:-1])))
    print('Writing images/%s' % (fname))
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
    if not os.path.exists('images'): os.system('mkdir -p images')
    if '/' in fname:
        os.system('mkdir -p images/%s' % ('/'.join(fname.split('/')[:-1])))
    print('Writing images/%s' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)
    plt.close()

def plot_spectra(source_fwd, source_adj):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    # Plot spectra
    ng = source_fwd.shape[1]
    xvals = np.array(range(ng))
    yvals = source_fwd[0, :]
    xvals = np.insert(xvals, [0, ng], [0, ng])
    yvals = np.insert(yvals, [0, ng], [0, 0 ])
    ax.step(xvals, yvals, where='post', label='Forward')
    for tally in source_adj.keys():
        yvals = source_adj[tally][0, :]
        yvals = np.insert(yvals, [0, ng], [0, 0])
        ax.step(xvals, yvals, where='post', label='Adjoint (%u)' % (tally))

    # Formatting
    ax.set_xlabel('Energy group')
    ax.set_ylabel('Magnitude')
    ax.set_title('Forward and Adjoint Source Energy Spectra')
    ax.grid()
    ax.legend()

    # Save figure
    if not os.path.exists('images'): os.makedirs('images')
    fname = 'spectra_lin.png'
    print('Writing images/%s' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)

    # Log scale
    ax.set_yscale('log')

    # Save figure
    fname = 'spectra_log.png'
    print('Writing images/%s' % (fname))
    plt.savefig('images/%s' % (fname), bbox_inches=bbox_inches)

    plt.close()

def get_mat_name_short(mat_name):
    if mat_name in mat_names_short.keys(): return mat_names_short[mat_name]
    return mat_name

mat_names_short = {
    'void'                               : 'Void'                 ,
    'Air (dry, near sea level)'          : 'Air'                  ,
    'Aluminum (Al)'                      : 'Al'                   ,
    'Aluminum Oxide (Al2O3)'             : 'Al2O3'                ,
    'Aluminum, alloy 6061-O'             : 'Al-6061O'             ,
    'Beryllium (Be)'                     : 'Be'                   ,
    'Beryllium Oxide (BeO)'              : 'BeO'                  ,
    'Bismuth (Bi)'                       : 'Bi'                   ,
    'Boron (B)'                          : 'B'                    ,
    'Cadmium (Cd)'                       : 'Cd'                   ,
    'Calcium Fluoride (CaF2)'            : 'CaF2'                 ,
    'Calcium Oxide (CaO)'                : 'CaO'                  ,
    'Carbon, Graphite (reactor grade)'   : 'Graphite'             ,
    'Concrete, M-1'                      : 'Concrete (M1)'        ,
    'Concrete, Regulatory Concrete'
    ' (developed for U.S. NRC)'          : 'Concrete (regulatory)',
    'Copper (Cu)'                        : 'Cu'                   ,
    'Earth, U.S. Average'                : 'Earth'                ,
    'Gadolinium (Gd)'                    : 'Gd'                   ,
    'Iron (Fe)'                          : 'Fe'                   ,
    'Lead (Pb)'                          : 'Pb'                   ,
    'Lithium (Li)'                       : 'Li'                   ,
    'Lithium Fluoride (LiF)'             : 'LiF'                  ,
    'Lithium Oxide (Li2O)'               : 'Li2O'                 ,
    'Magnesium (Mg)'                     : 'Mg'                   ,
    'Magnesium Oxide (MgO)'              : 'MgO'                  ,
    'Nickel (Ni)'                        : 'Ni'                   ,
    'Polyethylene, Borated'              : 'Borated HDPE'         ,
    'Polyethylene, Non-borated (C2H4)'   : 'HDPE'                 ,
    'Titanium (Ti)'                      : 'Ti'                   ,
    'Titanium Dioxide (TiO2)'            : 'TiO2'                 ,
    'Water, Heavy (D2O)'                 : 'D2O'                  ,
    'Water, Liquid (H2O)'                : 'H2O'                  ,
    'Zircaloy-4'                         : 'Zircaloy-4'           ,
    'Calcium (Ca)'                       : 'Ca'                   ,
    'Vanadium (V)'                       : 'V'                    ,
    'Nickel-60'                          : 'Ni-60'                ,
    'Magnesium Fluoride (MgF2)'          : 'MgF2'                 ,
    'Aluminum Fluoride (AlF3)'           : 'AlF3'                 ,
    'Titanium(III) Fluoride (TiF3)'      : 'TiF3'                 ,
    'Fluental'                           : 'Fluental'             ,
    '5% Borated Polyethylene (SWX-201)'  : '5% boron HDPE'        ,
    '30% Boron Polyethylene (SWX-210)'   : '30% boron HDPE'       ,
    '7.5% Lithium Polyethylene (SWX-215)': '7.5% lithium HDPE'    ,
    'Poly-Biz Gamma Shield (SWX-217)'    : 'Poly-Biz'             ,
    '70% AlF3, 30% Al'                   : '70% AlF3 30% Al'      }

if __name__ == '__main__': main()
