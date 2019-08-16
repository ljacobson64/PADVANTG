#!/usr/bin/python

import os
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from cycler import cycler

from draw_geom import draw_geometry

def main():
    data = read_all()
    plot_all(data)

def read_all():
    data = {}

    # Read pickles
    for fname in os.listdir('pickles'):
        if not fname.endswith('.npy'): continue
        if fname.startswith('angular'): continue
        data[fname[:-4]] = np.load('pickles/%s' % (fname))

    # Read binary dR
    with open('pickles/dR.bin', 'rb') as f: dR = np.fromfile(f, '<f8')
    dR_shape = (data['mat_names'].shape[0], data['material_map'].shape[0],
                data['material_map'].shape[1], data['material_map'].shape[2])
    data['dR'] = dR.reshape(dR_shape)

    return data

def plot_all(data):
    nm = data['mix_table'   ].shape[1]  # Number of pure materials
    nz = data['flux_fwd_int'].shape[0]  # Number of Z intervals
    
    g0 = data['mesh_g'][0]  # First energy group

    x_vals = data['mesh_x']  # Horizontal direction
    y_vals = data['mesh_y']  # Vertical direction
    
    hz = nz / 2  # Use the middle Z slice

    # Plot spectra
    plot_spectra(data['source_spectrum'], data['response_spectrum'])

    # Plot total source
    flux = np.sum(data['source'][hz, :, :, :], 2)
    title = 'Source (Total)'
    fname = 'source_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'lin')

    # Plot total response
    flux = np.sum(data['response'][hz, :, :, :], 2)
    title = 'Response (Total)'
    fname = 'response_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'lin')

    # Plot total forward flux
    flux = np.sum(data['flux_fwd'][hz, :, :, :], 2)
    title = 'Forward Flux (Total)'
    fname = 'flux_forward_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plot total adjoint flux
    flux = np.sum(data['flux_adj'][hz, :, :, :], 2)
    title = 'Adjoint Flux (Total)'
    fname = 'flux_adjoint_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plot total contributon
    flux = np.sum(data['contributon'][hz, :, :, :], 2)
    title = 'Contributon (Total)'
    fname = 'contributon_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plots for the fast group (2) and thermal group (26)
    for gx in [2, 26]:
        gf = gx - g0

        flux = data['flux_fwd'][hz, :, :, gf]
        title = 'Forward Flux (Group %u)' % (gx)
        fname = 'flux_forward_g%02u.png'  % (gx)
        plot_flux(flux, x_vals, y_vals, title, fname, 'log')

        flux = data['flux_adj'][hz, :, :, gf]
        title = 'Adjoint Flux (Group %u)' % (gx)
        fname = 'flux_adjoint_g%02u.png'  % (gx)
        plot_flux(flux, x_vals, y_vals, title, fname, 'log')

        flux = data['contributon'][hz, :, :, gf]
        title = 'Contributon (Group %u)' % (gx)
        fname = 'contributon_g%02u.png'  % (gx)
        plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plot dR for all materials
    for im in range(nm):
        flux = data['dR'][im, hz, :, :]
        title = r'$\delta R$ for Material %u (%s)' % (im, data['mat_names'][im])
        fname = 'dR_%02u.png' % (im)
        plot_flux(flux, x_vals, y_vals, title, fname, 'logplusminus')

def plot_flux(result, x_vals, y_vals, title, fname, fmt):
    plt.rcParams['font.size'       ] = 10.0
    plt.rcParams['font.family'     ] = 'serif'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['axes.prop_cycle' ] = cycler(color=['k'])

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8, 6)

    # Color scale
    if fmt == 'logplusminus':
        logvmax = 6
        loglinthresh = logvmax - 8
        vmax = 10**logvmax
        vmin = -vmax
        linthresh = 10**loglinthresh
        norm = colors.SymLogNorm(vmin=vmin, vmax=vmax,
                                 linthresh=linthresh, linscale=1.0)
        cmap='RdBu_r'
        ticks_neg = [-10**x for x in np.linspace(loglinthresh, logvmax, 9)]
        ticks_pos = [ 10**x for x in np.linspace(loglinthresh, logvmax, 9)]
        ticks = list(reversed(ticks_neg)) + [0.0] + ticks_pos
        tick_labels_neg = [r'$-10^{%u}$' % (x) for x in
                           np.arange(loglinthresh, logvmax + 1, 1)]
        tick_labels_pos = [r'$10^{%u}$'  % (x) for x in
                           np.arange(loglinthresh, logvmax + 1, 1)]
        tick_labels = (list(reversed(tick_labels_neg)) + [0] + tick_labels_pos)
    elif fmt == 'log':
        result_min_nonzero = np.min(result[result != 0])
        result_max_nonzero = np.max(result[result != 0])
        logvmin = np.floor(np.log10(result_min_nonzero))
        logvmax = np.ceil(np.log10(result_max_nonzero))
        vmin = 10**logvmin
        vmax = 10**logvmax
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        cmap = 'Spectral_r'
        ticks = [10**x for x in np.arange(logvmin, logvmax + 1, 1)]
        tick_labels = [r'$10^{%u}$' % (x) for x in
                       np.arange(logvmin, logvmax + 1, 1)]
    elif fmt == 'lin':
        vmin = 0.0
        vmax = np.max(result) * 1.1
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = 'Reds'

    # Plot data
    im = ax.pcolormesh(x_vals, y_vals, result, norm=norm, cmap=cmap)

    # Formatting
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
