#!/usr/bin/python

import os
import copy
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from cycler import cycler

from draw_geom import draw_geometry
from mcnp_colors import mcnp_colors as MC

def main():
    data = read_pickles()
    plot_all(data)

def read_pickles():
    data = {}
    for fname in os.listdir('pickles'):
        if not fname.endswith('.npy'): continue
        if fname.startswith('angular'): continue
        data[fname[:-4]] = np.load('pickles/%s' % (fname))
    return data

def plot_all(data):
    nm = data['mix_table'      ].shape[1]  # Number of pure materials
    nz = data['scalar_flux_fwd'].shape[0]  # Number of Z intervals
    
    g0 = data['mesh_g'][0]  # First energy group

    x_vals = data['mesh_x']  # Horizontal direction
    y_vals = data['mesh_y']  # Vertical direction
    
    hz = nz / 2  # Use the middle Z slice

    # Plot spectra
    plot_spectra(data['source_spectrum'], data['response_spectrum'])

    # Plot material map
    plot_data = data['material_map'][hz, :, :]
    title = 'Material Map'
    fname = 'material_map.png'
    plot_map(plot_data, x_vals, y_vals, title, fname, 'mats',
              mat_names=data['mat_names'])

    # Plot total source
    plot_data = np.sum(data['source'][hz, :, :, :], 2)
    title = 'Source (Total)'
    fname = 'source_total.png'
    plot_map(plot_data, x_vals, y_vals, title, fname, 'lin')

    # Plot total response
    plot_data = np.sum(data['response'][hz, :, :, :], 2)
    title = 'Response (Total)'
    fname = 'response_total.png'
    plot_map(plot_data, x_vals, y_vals, title, fname, 'lin')

    # Plot total forward flux
    plot_data = np.sum(data['scalar_flux_fwd'][hz, :, :, :], 2)
    title = 'Forward Flux (Total)'
    fname = 'scalar_flux_fwd_total.png'
    plot_map(plot_data, x_vals, y_vals, title, fname, 'log')

    # Plot total adjoint flux
    plot_data = np.sum(data['scalar_flux_adj'][hz, :, :, :], 2)
    title = 'Adjoint Flux (Total)'
    fname = 'scalar_flux_adj_total.png'
    plot_map(plot_data, x_vals, y_vals, title, fname, 'log')

    # Plot total contributon flux
    plot_data = np.sum(data['scalar_flux_con'][hz, :, :, :], 2)
    title = 'Contributon Flux (Total)'
    fname = 'scalar_flux_con_total.png'
    plot_map(plot_data, x_vals, y_vals, title, fname, 'log')

    # Plot total forward current
    plot_data = np.sum(data['current_fwd'][hz, :, :, :, :2], 2)
    title = 'Forward Current (Total)'
    fname = 'current_fwd_total.png'
    plot_quiver(plot_data, x_vals, y_vals, title, fname)

    # Plot total adjoint current
    plot_data = np.sum(data['current_adj'][hz, :, :, :, :2], 2)
    title = 'Adjoint Current (Total)'
    fname = 'current_adj_total.png'
    plot_quiver(plot_data, x_vals, y_vals, title, fname)

    # Plot total contributon current
    plot_data = np.sum(data['current_con'][hz, :, :, :, :2], 2)
    title = 'Contributon Current (Total)'
    fname = 'current_con_total.png'
    plot_quiver(plot_data, x_vals, y_vals, title, fname)

    # Plots for the fast group (2) and thermal group (26)
    for igx in [2, 26]:
        igf = igx - g0

        # Forward flux
        plot_data = data['scalar_flux_fwd'][hz, :, :, igf]
        title = 'Forward Flux (Group %u)' % (igx)
        fname = 'scalar_flux_fwd_g%02u.png'  % (igx)
        plot_map(plot_data, x_vals, y_vals, title, fname, 'log')

        # Adjoint flux
        plot_data = data['scalar_flux_adj'][hz, :, :, igf]
        title = 'Adjoint Flux (Group %u)' % (igx)
        fname = 'scalar_flux_adj_g%02u.png'  % (igx)
        plot_map(plot_data, x_vals, y_vals, title, fname, 'log')

        # Contributon flux
        plot_data = data['scalar_flux_con'][hz, :, :, igf]
        title = 'Contributon Flux (Group %u)' % (igx)
        fname = 'scalar_flux_con_g%02u.png'  % (igx)
        plot_map(plot_data, x_vals, y_vals, title, fname, 'log')

        # Forward current
        plot_data = data['current_fwd'][hz, :, :, igf, :2]
        title = 'Forward Current (Group %u)' % (igx)
        fname = 'current_fwd_g%02u.png'  % (igx)
        plot_quiver(plot_data, x_vals, y_vals, title, fname)

        # Adjoint current
        plot_data = data['current_adj'][hz, :, :, igf, :2]
        title = 'Adjoint Current (Group %u)' % (igx)
        fname = 'current_adj_g%02u.png'  % (igx)
        plot_quiver(plot_data, x_vals, y_vals, title, fname)

        # Contributon current
        plot_data = data['current_con'][hz, :, :, igf, :2]
        title = 'Contributon Current (Group %u)' % (igx)
        fname = 'current_con_g%02u.png'  % (igx)
        plot_quiver(plot_data, x_vals, y_vals, title, fname)

    # Plot dR for all materials
    for im in range(nm):
        plot_data = data['dR_angular'][im, hz, :, :]
        title = (r'$\delta R$ (Angular) for Material %u (%s)' %
                 (im, data['mat_names'][im]))
        fname = 'dR_angular_%02u.png' % (im)
        plot_map(plot_data, x_vals, y_vals, title, fname, 'logplusminus')
    for im in range(nm):
        plot_data = data['dR_scalar'][im, hz, :, :]
        title = (r'$\delta R$ (Scalar) for Material %u (%s)' %
                 (im, data['mat_names'][im]))
        fname = 'dR_scalar_%02u.png' % (im)
        plot_map(plot_data, x_vals, y_vals, title, fname, 'logplusminus')

def plot_map(plot_data, x_vals, y_vals, title, fname, fmt, mat_names=None):
    plt.rcParams['font.size'       ] = 10.0
    plt.rcParams['font.family'     ] = 'serif'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['axes.prop_cycle' ] = cycler(color=['k'])

    # Color scale
    if fmt == 'logplusminus':
        logvmax = 8
        span = 8
        loglinthresh = logvmax - span
        vmax = 10**logvmax
        vmin = -vmax
        linthresh = 10**loglinthresh
        norm = colors.SymLogNorm(vmin=vmin, vmax=vmax,
                                 linthresh=linthresh, linscale=1.0)
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
    elif fmt == 'log':
        result_min_nonzero = np.min(plot_data[plot_data != 0])
        result_max_nonzero = np.max(plot_data[plot_data != 0])
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
        vmax = np.max(plot_data) * 1.1
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = 'Reds'
    elif fmt == 'mats':
        with open('com_1', 'r') as reader: lines = reader.readlines()
        CL = [MC['white']]
        for line in lines:
            if not line.startswith('shade'): continue
            tokens = line.split()
            CL.append(MC[tokens[2]])
        vmin = 0
        vmax = len(CL)
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(CL)
        ticks = np.linspace(0.5, vmax - 0.5, len(CL))
        tick_labels = list(mat_names)

    fig, ax = plt.subplots(1, 1)
    if 'mats' in fmt:
        fig.set_size_inches(10, 6)
        plt.subplots_adjust(left=-0.025, right=0.75)
    else:
        fig.set_size_inches(8, 6)
        plt.subplots_adjust(left=0.125, right=0.9)

    # Plot data
    im = ax.pcolormesh(x_vals, y_vals, plot_data, norm=norm, cmap=cmap)

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
    #plt.tight_layout()
    plt.savefig('images/%s' % (fname))
    plt.close()

def plot_quiver(plot_data, x_vals, y_vals, title, fname):
    plt.rcParams['font.size'       ] = 10.0
    plt.rcParams['font.family'     ] = 'serif'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['axes.prop_cycle' ] = cycler(color=['k'])

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
    #plt.tight_layout()
    plt.savefig('images/%s' % (fname))
    plt.close()

def plot_spectra(source, response):
    plt.rcParams['font.size'       ] = 10.0
    plt.rcParams['font.family'     ] = 'serif'
    plt.rcParams['mathtext.default'] = 'regular'

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    # Plot spectra
    ng = len(source)
    ebins = range(ng)
    plt.bar(ebins, source  , label='Source'  )
    plt.bar(ebins, response, label='Response')

    # Formatting
    ax.set_xlim([-0.6, ng - 0.4])
    ax.set_xlabel('Energy group')
    ax.set_ylabel('Magnitude')
    ax.set_title('Source and Response Energy Spectra')
    ax.legend()

    # Save figure
    if not os.path.exists('images'): os.makedirs('images')
    fname = 'spectra_lin.png'
    print('Writing images/%s to file' % (fname))
    #plt.tight_layout()
    plt.savefig('images/%s' % (fname))

    # Log scale
    ax.set_yscale('log')

    # Save figure
    fname = 'spectra_log.png'
    print('Writing images/%s to file' % (fname))
    #plt.tight_layout()
    plt.savefig('images/%s' % (fname))

    plt.close()

if __name__ == '__main__': main()
