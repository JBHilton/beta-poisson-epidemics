import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import arange, array, ceil, floor, isnan, linspace, log10, where
from pickle import load
from scipy import stats
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)

formats = ['.png', '.svg', '.eps']

data_name_list = [
    'plague_data',
    'monkeypox_data',
    'fasina_ebola_data',
    'fay_ebola_data',
    'cdc_sars_data',
    'cowling_mers_data',
    'mers_data',
    'noro_data'
    ]
data_set_list = [
    plague_data,
    monkeypox_data,
    fasina_ebola_data,
    fay_ebola_data,
    cdc_sars_data,
    cowling_mers_data,
    mers_data,
    noro_data
    ]

mle_list = []
for i, data_name in enumerate(data_name_list):
    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname, 'rb') as f:
        (mle_dict,
        var_dict,
        ci_dict,
        var_ci_dict,
        llh_dict) = load(f)
    mle_list.append(mle_dict)

lmbd_vals_list = []
phi_vals_list = []
nu_vals_list = []
lmbd_curve_list = []
phi_curve_list = []
nu_curve_list = []
lmbd_grid_list = []
phi_grid_list = []
nu_grid_list = []
for i, data_name in enumerate(data_name_list):
    fname = 'outputs/sensitivity_analyses/'+data_name+'_results.pkl'
    with open(fname, 'rb') as f:
        (lmbd_vals,
        phi_vals,
        nu_vals,
        lmbd_curve,
        phi_curve,
        nu_curve,
        lmbd_grid,
        phi_grid,
        nu_grid) = load(f)
    lmbd_vals_list.append(lmbd_vals)
    phi_vals_list.append(phi_vals)
    nu_vals_list.append(nu_vals)
    lmbd_curve_list.append(lmbd_curve)
    phi_curve_list.append(phi_curve)
    nu_curve_list.append(nu_curve)
    lmbd_grid_list.append(lmbd_grid)
    phi_grid_list.append(phi_grid)
    nu_grid_list.append(nu_grid)

def calculate_label_pos(ax):
    xlims = ax.get_xlim()
    xrange = xlims[1]-xlims[0]
    ylims = ax.get_ylim()
    yrange = ylims[1]-ylims[0]

    return [xlims[0]-0.3*xrange, ylims[1]+0.05*yrange]

for i in range(len(data_name_list)):

    true_lmbd = mle_list[i]['beta-Poisson'][0]
    true_phi = mle_list[i]['beta-Poisson'][1]
    true_nu = mle_list[i]['beta-Poisson'][2]

    no_ticks = 5

    phi_grid_nonzero = where(~isnan(phi_grid_list[i]))
    nu_grid_nonzero = where(~isnan(nu_grid_list[i]))



    fig, (lmbd_ax, phi_ax, nu_ax) = plt.subplots(1,
                                                 3,
                                                 constrained_layout=True)

    lmbd_ax.plot(lmbd_vals_list[i], lmbd_curve_list[i], 'k')
    lmbd_ax.plot([true_lmbd, true_lmbd], lmbd_ax.get_ylim(), '--b')
    lmbd_ax.set_aspect(1.0/lmbd_ax.get_data_ratio())
    lmbd_ax.set_xlabel('$\\lambda$')
    lmbd_ax.set_ylabel('Log likelihood')
    label_pos = calculate_label_pos(lmbd_ax)
    lmbd_ax.text(label_pos[0], label_pos[1], 'a)',
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    phi_ax.plot(phi_vals_list[i], phi_curve_list[i], 'k')
    phi_ax.plot([true_phi, true_phi], phi_ax.get_ylim(), '--b')
    phi_ax.set_aspect(1.0/phi_ax.get_data_ratio())
    phi_ax.set_xlabel('$\\Phi$')
    phi_ax.set_ylabel('Log likelihood')
    label_pos = calculate_label_pos(phi_ax)
    phi_ax.text(label_pos[0], label_pos[1], 'b)',
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    nu_ax.plot(nu_vals_list[i], nu_curve_list[i], 'k')
    nu_ax.plot([true_nu, true_nu], nu_ax.get_ylim(), '--b')
    nu_ax.set_aspect(1.0/nu_ax.get_data_ratio())
    nu_ax.set_xlabel('$\\nu$')
    nu_ax.set_ylabel('Log likelihood')
    label_pos = calculate_label_pos(nu_ax)
    nu_ax.text(label_pos[0], label_pos[1], 'c)',
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    for fmt in formats:
        fig.savefig('plots/'+data_name_list[i]+'_sensitivity_curves'+fmt, bbox_inches='tight')
    plt.close()

    fig, (lmbd_ax, phi_ax, nu_ax) = plt.subplots(1,
                                                 3,
                                                 constrained_layout=True,
                                                 figsize=(9, 3))

    lmbd_window = lmbd_vals_list[i][-1] - lmbd_vals_list[i][0]
    phi_window = phi_vals_list[i][-1] - phi_vals_list[i][0]
    nu_window = nu_vals_list[i][-1] - nu_vals_list[i][0]

    axim=lmbd_ax.imshow(lmbd_grid_list[i],
           origin='lower',
           extent=(phi_vals_list[i][0],
                   phi_vals_list[i][-1],
                   nu_vals_list[i][0],
                   nu_vals_list[i][-1]),
           aspect='auto')
    lmbd_ax.plot(true_phi, true_nu, 'kx', markersize=5)
    lmbd_ax.text(true_phi+0.01*phi_window, true_nu+0.02*nu_window, '($\hat{\Phi}$, $\hat{\\nu}$)', fontsize=9)
    lmbd_ax.set_xlabel('$\\Phi$')
    lmbd_ax.set_ylabel('$\\nu$')
    lmbd_ax.text(-0.25 * phi_window, 1.1 * nu_window, 'a)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
    divider = make_axes_locatable(lmbd_ax)
    cax = divider.append_axes("top", size="5%", pad=.25)
    cbar = plt.colorbar(axim,
                    label="Log likelihood",
                    cax=cax,
                    orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='x', rotation=45)
    lmbd_ax.spines['top'].set_visible(False)
    lmbd_ax.spines['right'].set_visible(False)
    cbar.outline.set_visible(False)

    axim=phi_ax.imshow(phi_grid_list[i],
           origin='lower',
           extent=(lmbd_vals_list[i][0],
                   lmbd_vals_list[i][-1],
                   nu_vals_list[i][0],
                   nu_vals_list[i][-1]),
           aspect='auto')
    phi_ax.plot(true_lmbd, true_nu, 'kx', markersize=5)
    phi_ax.text(true_lmbd+0.01*lmbd_window, true_nu+0.02*nu_window, '($\hat{\lambda}$, $\hat{\\nu}$)', fontsize=9)
    phi_ax.set_xlabel('$\\lambda$')
    phi_ax.set_ylabel('$\\nu$')
    phi_ax.text(-0.25 * lmbd_window, 1.1 * nu_window, 'b)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
    divider = make_axes_locatable(phi_ax)
    cax = divider.append_axes("top", size="5%", pad=.25)
    cbar = plt.colorbar(axim,
                    label="Log likelihood",
                    cax=cax,
                    orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='x', rotation=45)
    phi_ax.spines['top'].set_visible(False)
    phi_ax.spines['right'].set_visible(False)
    cbar.outline.set_visible(False)

    axim=nu_ax.imshow(nu_grid_list[i],
           origin='lower',
           extent=(lmbd_vals_list[i][0],
                   lmbd_vals_list[i][-1],
                   phi_vals_list[i][0],
                   phi_vals_list[i][-1]),
           aspect='auto')
    nu_ax.plot(true_lmbd, true_phi, 'kx', markersize=5)
    nu_ax.text(true_lmbd+0.01*lmbd_window, true_phi+0.02*phi_window, '($\hat{\lambda}$, $\hat{\Phi}$)', fontsize=9)
    nu_ax.set_xlabel('$\\lambda$')
    nu_ax.set_ylabel('$\\Phi$')
    nu_ax.text(-0.25 * lmbd_window, 1.1 * phi_window, 'c)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
    divider = make_axes_locatable(nu_ax)
    cax = divider.append_axes("top", size="5%", pad=.25)
    cbar = plt.colorbar(axim,
                    label="Log likelihood",
                    cax=cax,
                    orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='x', rotation=45)
    nu_ax.spines['top'].set_visible(False)
    nu_ax.spines['right'].set_visible(False)
    cbar.outline.set_visible(False)

    for fmt in formats:
        fig.savefig('plots/'+data_name_list[i]+'_sensitivity_grids'+fmt, bbox_inches='tight')
    plt.close()

    lmbd_grid_scale = 10**floor(log10(-lmbd_grid_list[i].mean()))
    phi_grid_scale = 10**floor(log10(-phi_grid_list[i][phi_grid_nonzero].mean()))
    nu_grid_scale = 10**floor(log10(-nu_grid_list[i][nu_grid_nonzero].mean()))

    lmbd_grid_min = lmbd_grid_scale * \
                    floor(lmbd_grid_list[i].min() / lmbd_grid_scale)
    lmbd_grid_max = lmbd_grid_scale * \
                    ceil(lmbd_grid_list[i].max() / lmbd_grid_scale)
    lmbd_grid_tick = arange(lmbd_grid_min, lmbd_grid_max+lmbd_grid_scale, lmbd_grid_scale)
    if len(lmbd_grid_tick)>10:
        no_10s = int(len(lmbd_grid_tick)/10)
        lmbd_grid_tick = lmbd_grid_tick[::2**no_10s]
    elif len(lmbd_grid_tick)<5:
        lmbd_grid_tick = arange(lmbd_grid_min, lmbd_grid_max, 0.5 * lmbd_grid_scale)

    phi_grid_min = phi_grid_scale * \
                    floor(phi_grid_list[i][phi_grid_nonzero].min() / phi_grid_scale)
    phi_grid_max = phi_grid_scale * \
                    ceil(phi_grid_list[i][phi_grid_nonzero].max() / phi_grid_scale)
    phi_grid_tick = arange(phi_grid_min, phi_grid_max+phi_grid_scale, phi_grid_scale)
    if len(phi_grid_tick)>10:
        no_10s = int(len(phi_grid_tick)/10)
        phi_grid_tick = phi_grid_tick[::2**no_10s]
    elif len(phi_grid_tick)<5:
        phi_grid_tick = arange(phi_grid_min, phi_grid_max, 0.5 * phi_grid_scale)

    nu_grid_min = nu_grid_scale * \
                    floor(nu_grid_list[i][nu_grid_nonzero].min() / nu_grid_scale)
    nu_grid_max = nu_grid_scale * \
                    ceil(nu_grid_list[i][nu_grid_nonzero].max() / nu_grid_scale)
    nu_grid_tick = arange(nu_grid_min, nu_grid_max+nu_grid_scale, nu_grid_scale)
    if len(nu_grid_tick)>10:
        no_10s = int(len(lmbd_grid_tick)/10)
        nu_grid_tick = nu_grid_tick[::2**no_10s]
    elif len(nu_grid_tick)<5:
        nu_grid_tick = arange(nu_grid_min, nu_grid_max, 0.5 * nu_grid_scale)

    fig, (lmbd_ax, phi_ax, nu_ax) = plt.subplots(1,
                                                 3,
                                                 constrained_layout=True,
                                                 figsize=(9, 3))

    lmbd_lvs = lmbd_grid_tick
    phi_lvs = phi_grid_tick
    nu_lvs = nu_grid_tick

    axim=lmbd_ax.contour(lmbd_grid_list[i],
           origin='lower',
           colors='k',
           levels=lmbd_lvs,
           extent=(phi_vals_list[i][0],
                   phi_vals_list[i][-1],
                   nu_vals_list[i][0],
                   nu_vals_list[i][-1]),
           vmin=lmbd_grid_min,
           vmax=lmbd_grid_max)
    lmbd_ax.clabel(axim, fontsize=9, inline=1, fmt='%1.1f')
    lmbd_ax.set_xlabel('$\\Phi$')
    lmbd_ax.set_ylabel('$\\nu$')
    lmbd_ax.text(-0.25 * phi_window, 1.1 * nu_window, 'a)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    axim=phi_ax.contour(phi_grid_list[i],
           origin='lower',
           colors='k',
           levels=phi_lvs,
           extent=(lmbd_vals_list[i][0],
                   lmbd_vals_list[i][-1],
                   nu_vals_list[i][0],
                   nu_vals_list[i][-1]),
           vmin=phi_grid_min,
           vmax=phi_grid_max)
    phi_ax.clabel(axim, fontsize=9, inline=1, fmt='%1.1f')
    phi_ax.set_xlabel('$\\lambda$')
    phi_ax.set_ylabel('$\\nu$')
    phi_ax.text(-0.25 * lmbd_window, 1.1 * nu_window, 'b)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    axim=nu_ax.contour(nu_grid_list[i],
           origin='lower',
           colors='k',
           levels=nu_lvs,
           extent=(lmbd_vals_list[i][0],
                   lmbd_vals_list[i][-1],
                   phi_vals_list[i][0],
                   phi_vals_list[i][-1]),
           vmin=nu_grid_min,
           vmax=nu_grid_max)
    nu_ax.clabel(axim, fontsize=9, inline=1, fmt='%1.1f')
    nu_ax.set_xlabel('$\\lambda$')
    nu_ax.set_ylabel('$\\Phi$')
    nu_ax.text(-0.25 * lmbd_window, 1.1 * phi_window, 'c)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    for fmt in formats:
        fig.savefig('plots/'+data_name_list[i]+'_sensitivity_contours'+fmt, bbox_inches='tight')
    plt.close()
