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
    ng_tot   = data['sigma_t'  ].shape[1]
    ng_flux  = data['flux_fwd' ].shape[0]
    nx       = data['flux_fwd' ].shape[1]
    ny       = data['flux_fwd' ].shape[2]
    nz       = data['flux_fwd' ].shape[3]
    num_mats = data['mix_table'].shape[1]

    g0 = data['mesh_g'][0]

    x_vals = data['mesh_x']
    y_vals = data['mesh_y']
    hz = nz / 2

    # Plot spectra
    plot_spectra(data['source_spectrum'], data['response_spectrum'])

    # Plot total source
    flux = np.flipud(np.rot90(np.sum(data['source'][:, :, :, hz], 0)))
    title = 'Source (Total)'
    fname = 'source_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'lin')

    # Plot total response
    flux = np.flipud(np.rot90(np.sum(data['response'][:, :, :, hz], 0)))
    title = 'Response (Total)'
    fname = 'response_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'lin')

    # Plot total forward flux
    flux = np.flipud(np.rot90(np.sum(data['flux_fwd'][:, :, :, hz], 0)))
    title = 'Forward Flux (Total)'
    fname = 'flux_forward_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plot total adjoint flux
    flux = np.flipud(np.rot90(np.sum(data['flux_adj'][:, :, :, hz], 0)))
    title = 'Adjoint Flux (Total)'
    fname = 'flux_adjoint_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plot total contributon
    flux = np.flipud(np.rot90(np.sum(data['contributon'][:, :, :, hz], 0)))
    title = 'Contributon (Total)'
    fname = 'contributon_total.png'
    plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plots for the fast group (2) and thermal group (26)
    for g in [2, 26]:
        flux = np.flipud(np.rot90(data['flux_fwd'][g - g0, :, :, hz]))
        title = 'Forward Flux (Group %u)' % (g)
        fname = 'flux_forward_g%02u.png' % (g)
        plot_flux(flux, x_vals, y_vals, title, fname, 'log')

        flux = np.flipud(np.rot90(data['flux_adj'][g - g0, :, :, hz]))
        title = 'Adjoint Flux (Group %u)' % (g)
        fname = 'flux_adjoint_g%02u.png' % (g)
        plot_flux(flux, x_vals, y_vals, title, fname, 'log')

        flux = np.flipud(np.rot90(data['contributon'][g - g0, :, :, hz]))
        title = 'Contributon (Group %u)' % (g)
        fname = 'contributon_g%02u.png' % (g)
        plot_flux(flux, x_vals, y_vals, title, fname, 'log')

    # Plot dR for all materials
    for i in range(num_mats):
        flux = np.flipud(np.rot90(data['dR'][i, :, :, hz]))
        title = r'$\delta R$ for Material %u (%s)' % (i, data['mat_names'][i])
        fname = 'dR_%02u.png' % (i)
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
