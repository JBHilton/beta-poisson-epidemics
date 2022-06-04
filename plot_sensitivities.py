import matplotlib.pyplot as plt
from numpy import arange
from pickle import load
from scipy import stats
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)

data_name_list = [
    'plague_data'
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

lmbd_vals = arange(.1, 5., .1)
phi_vals = arange(.1, 10., .1)
N_inv_vals = arange(0., 1., .1)

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

lmbd_curve_list = []
phi_curve_list = []
N_inv_curve_list = []
lmbd_grid_list = []
phi_grid_list = []
N_inv_grid_list = []
for i, data_name in enumerate(data_name_list):
    fname = 'outputs/sensitivity_analyses/'+data_name+'_results.pkl'
    with open(fname, 'rb') as f:
        (lmbd_curve,
        phi_curve,
        N_inv_curve,
        lmbd_grid,
        phi_grid,
        N_inv_grid) = load(f)
    lmbd_curve_list.append(lmbd_curve)
    phi_curve_list.append(phi_curve)
    N_inv_curve_list.append(N_inv_curve)
    lmbd_grid_list.append(lmbd_grid)
    phi_grid_list.append(phi_grid)
    N_inv_grid_list.append(N_inv_grid)

for i in range(len(data_name_list)):

    fig, (lmbd_ax, phi_ax, N_inv_ax) =plt.subplots(1, 3)

    print(len(lmbd_curve_list[i]))

    true_lmbd = mle_list[i]['beta-Poisson'][0]
    true_phi = mle_list[i]['beta-Poisson'][1]
    true_N_inv = mle_list[i]['beta-Poisson'][2]

    lmbd_ax.plot(lmbd_vals, lmbd_curve_list[i], 'k', linewidth=6)
    lmbd_ax.plot([true_lmbd, true_lmbd], lmbd_ax.get_ylim(), '--b', linewidth=6)
    lmbd_ax.set_xlabel('$\\lambda$', fontsize=40)
    lmbd_ax.set_ylabel('Log likelihood', fontsize=40)
    lmbd_ax.set_xticks(fontsize=40)
    lmbd_ax.set_yticks(fontsize=40)

    phi_ax.plot(phi_vals, phi_curve_list[i], 'k', linewidth=6)
    phi_ax.plot([true_phi, true_phi], phi_ax.get_ylim(), '--b', linewidth=6)
    phi_ax.set_xlabel('$\\Phi$', fontsize=40)
    phi_ax.set_ylabel('Log likelihood', fontsize=40)
    phi_ax.set_xticks(fontsize=40)
    phi_ax.set_yticks(fontsize=40)

    N_inv_ax.plot(N_inv_vals, N_inv_curve_list[i], 'k', linewidth=6)
    N_inv_ax.plot([true_N_inv, true_N_inv], N_inv_ax.get_ylim(), '--b', linewidth=6)
    N_inv_ax.set_xlabel('$\\nu$', fontsize=40)
    N_inv_ax.set_ylabel('Log likelihood', fontsize=40)
    N_inv_ax.set_xticks(fontsize=30)
    N_inv_ax.set_yticks(fontsize=40)

    fig.savefig(data_name_list[i]+'sensitivity_curves.png', format='png', bbox_inches='tight')
