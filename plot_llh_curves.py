"This script generates plots of the log likelihood curves around the MLEs of the beta-Poisson "
"parameters for each dataset."

from os import mkdir
from os.path import isdir
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import arange, array, ceil, floor, isnan, linspace, log10, where
from pickle import load
from scipy import stats
from datasets import (plague_data, mpox_data, nigeria_ebola_data,
    guinea_ebola_data, singapore_sars_data, sk_mers_data, sa_mers_data, noro_data)

if isdir('outputs/sensitivity_analyses'):
    fname_root = 'outputs/'
else:
    fname_root = 'reference-outputs/'

if isdir('plots') is False:
    mkdir('plots')

plt.rcParams.update({'font.size': 16})

formats = ['.png', '.svg', '.eps']

data_key_list = [
    'plague_data',
    'mpox_data',
    'nigeria_ebola_data',
    'guinea_ebola_data',
    'singapore_sars_data',
    'sk_mers_data',
    'sa_mers_data',
    'noro_data'
    ]
data_name_list = [
    'Plague',
    'Mpox',
    'Ebola, Nigeria 2014',
    'Ebola, Guinea 2014',
    'SARS, Singapore 2003',
    'MERS, South Korea 2015',
    'MERS, Saudi Arabia 2015',
    'Norovirus, Netherlands 2012'
    ]
data_set_list = [
    plague_data,
    mpox_data,
    nigeria_ebola_data,
    guinea_ebola_data,
    singapore_sars_data,
    sk_mers_data,
    sa_mers_data,
    noro_data
    ]
figlabels = ['a)',
    'b)',
    'c)',
    'd)',
    'e)',
    'f)',
    'g)',
    'h)']

mle_list = []
ci_list = []
superspread_list = []
superspread_ci_list = []
p0_list = []
p0_ci_list = []
od_list = []
od_ci_list = []
for i, data_key in enumerate(data_key_list):
    fname = fname_root+'mles/'+data_key+'_results.pkl'
    with open(fname,'rb') as f:
        (mle_dict,
            var_dict,
            od_dict,
            superspread_dict,
            p0_dict,
            ci_dict,
            mean_ci_dict,
            var_ci_dict,
            od_ci_dict,
            superspread_ci_dict,
            p0_ci_dict,
            llh_dict) = load(f)
    mle_list.append(mle_dict)
    ci_list.append(ci_dict)
    superspread_list.append(superspread_dict)
    superspread_ci_list.append(superspread_ci_dict)
    p0_list.append(p0_dict)
    p0_ci_list.append(p0_ci_dict)
    od_list.append(od_dict)
    od_ci_list.append(od_ci_dict)

lmbd_vals_list = []
phi_vals_list = []
nu_vals_list = []
lmbd_curve_list = []
phi_curve_list = []
nu_curve_list = []
lmbd_grid_list = []
phi_grid_list = []
nu_grid_list = []
for i, data_key in enumerate(data_key_list):
    fname = fname_root+'sensitivity_analyses/'+data_key+'_results.pkl'
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

fig, axes = plt.subplots(8, 3, figsize=(20, 35))
fig.tight_layout()
plt.subplots_adjust(hspace=0.3)
axes = axes.flatten()

for i in range(len(data_name_list)):

    true_lmbd = mle_list[i]['beta-Poisson'][0]
    true_phi = mle_list[i]['beta-Poisson'][1]
    true_nu = mle_list[i]['beta-Poisson'][2]

    no_ticks = 5

    phi_grid_nonzero = where(~isnan(phi_grid_list[i]))
    nu_grid_nonzero = where(~isnan(nu_grid_list[i]))

    axes[i*3].plot(lmbd_vals_list[i], lmbd_curve_list[i], 'k')
    axes[i*3].plot([true_lmbd, true_lmbd], axes[i*3].get_ylim(), '--b')
    axes[i*3].set_aspect(1.0/axes[i*3].get_data_ratio())
    axes[i*3].set_xlabel('$\\lambda$')
    axes[i*3].set_ylabel('Log likelihood')
    label_pos = calculate_label_pos(axes[i*3])
    axes[i*3].text(label_pos[0], label_pos[1], figlabels[i]+"(i)",
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    axes[i*3+1].plot(phi_vals_list[i], phi_curve_list[i], 'k')
    axes[i*3+1].plot([true_phi, true_phi], axes[i*3+1].get_ylim(), '--b')
    axes[i*3+1].set_aspect(1.0/axes[i*3+1].get_data_ratio())
    axes[i*3+1].set_xlabel('$\\Phi$')
    axes[i*3+1].set_ylabel('Log likelihood')
    label_pos = calculate_label_pos(axes[i*3+1])
    axes[i*3+1].text(label_pos[0], label_pos[1], figlabels[i]+"(ii)",
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

    axes[i*3+2].plot(nu_vals_list[i], nu_curve_list[i], 'k', label="Log likelihood")
    axes[i*3+2].plot([true_nu, true_nu], axes[i*3+2].get_ylim(), '--b', label="Location of MLE")
    axes[i*3+2].set_aspect(1.0/axes[i*3+2].get_data_ratio())
    axes[i*3+2].set_xlabel('$\\nu$')
    axes[i*3+2].set_ylabel('Log likelihood')
    label_pos = calculate_label_pos(axes[i*3+2])
    axes[i*3+2].text(label_pos[0], label_pos[1], figlabels[i]+"(iii)",
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

axes[2].legend(loc='center left', bbox_to_anchor=(1.1, 0.5))

for fmt in formats:
    fig.savefig('plots/llh_curves'+fmt, bbox_inches='tight')
plt.close()