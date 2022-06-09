import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import arange, array, ceil, floor, isnan, where
from pickle import load
from scipy import stats
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)

data_name_list = [
    'plague_data',
    'monkeypox_data'
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

    lmbd_grid_min = 5 * floor(0.2 * lmbd_grid_list[i].min())
    lmbd_grid_max = 5 * ceil(0.2 * lmbd_grid_list[i].max())
    lmbd_grid_tick = lmbd_grid_min + 1/(no_ticks-1) * \
                    (lmbd_grid_max - lmbd_grid_min) * arange(no_ticks)
    phi_grid_min = 5 * floor(0.2 * phi_grid_list[i][phi_grid_nonzero].min())
    phi_grid_max = 5 * ceil(0.2 * phi_grid_list[i][phi_grid_nonzero].max())
    phi_grid_tick = phi_grid_min + 1/(no_ticks-1) * \
                    (phi_grid_max - phi_grid_min) * arange(no_ticks)
    nu_grid_min = 5 * floor(0.2 * nu_grid_list[i][nu_grid_nonzero].min())
    nu_grid_max = 5 * ceil(0.2 * nu_grid_list[i][nu_grid_nonzero].max())
    nu_grid_tick = nu_grid_min + 1/(no_ticks-1) * \
                    (nu_grid_max - nu_grid_min) * arange(no_ticks)

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

    fig.savefig(data_name_list[i]+'_sensitivity_curves.png', bbox_inches='tight')

    fig, (lmbd_ax, phi_ax, nu_ax) = plt.subplots(1,
                                                 3,
                                                 constrained_layout=True)

    axim=lmbd_ax.imshow(lmbd_grid_list[i],
           origin='lower',
           extent=(phi_vals_list[i][0],
                   phi_vals_list[i][-1],
                   nu_vals_list[i][0],
                   nu_vals_list[i][-1]),
           vmin=lmbd_grid_min,
           vmax=lmbd_grid_max,
           aspect='auto')
    lmbd_ax.set_xlabel('$\\Phi$')
    lmbd_ax.set_ylabel('$\\nu$')
    lmbd_ax.text(-0.5, 1.1, 'a)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
    divider = make_axes_locatable(lmbd_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(axim,
                    label="Log likelihood",
                    cax=cax,
                    ticks=lmbd_grid_tick)
    lmbd_ax.spines['top'].set_visible(False)
    lmbd_ax.spines['right'].set_visible(False)
    cbar.outline.set_visible(False)

    axim=phi_ax.imshow(phi_grid_list[i],
           origin='lower',
           extent=(lmbd_vals_list[i][0],
                   lmbd_vals_list[i][-1],
                   nu_vals_list[i][0],
                   nu_vals_list[i][-1]),
           vmin=phi_grid_min,
           vmax=phi_grid_max,
           aspect='auto')
    phi_ax.set_xlabel('$\\lambda$')
    phi_ax.set_ylabel('$\\nu$')
    phi_ax.text(-0.5, 1.1, 'b)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
    divider = make_axes_locatable(phi_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(axim,
                    label="Log likelihood",
                    cax=cax,
                    ticks=phi_grid_tick)
    phi_ax.spines['top'].set_visible(False)
    phi_ax.spines['right'].set_visible(False)
    cbar.outline.set_visible(False)

    axim=nu_ax.imshow(nu_grid_list[i],
           origin='lower',
           extent=(lmbd_vals_list[i][0],
                   lmbd_vals_list[i][-1],
                   phi_vals_list[i][0],
                   phi_vals_list[i][-1]),
           vmin=nu_grid_min,
           vmax=nu_grid_max,
           aspect='auto')
    nu_ax.set_xlabel('$\\lambda$')
    nu_ax.set_ylabel('$\\Phi$')
    nu_ax.text(-0.5, 1.1, 'c)',
                fontsize='medium', verticalalignment='top', fontfamily='serif',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
    divider = make_axes_locatable(nu_ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(axim,
                    label="Log likelihood",
                    cax=cax,
                    ticks=nu_grid_tick)
    nu_ax.spines['top'].set_visible(False)
    nu_ax.spines['right'].set_visible(False)
    cbar.outline.set_visible(False)

    fig.savefig(data_name_list[i]+'_sensitivity_grids.png', bbox_inches='tight')
