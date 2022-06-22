'''This plots the output of the plague analysis.'''

import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)
from functions import beta_poisson_pmf, zip_pmf

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
figlabels = ['a)',
    'b)',
    'c)',
    'd)',
    'e)',
    'f)',
    'g)',
    'h)']

mle_list = []
for i, data_name in enumerate(data_name_list):
    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname,'rb') as f:
        (mle_dict,
        var_dict,
        ci_dict,
        var_ci_dict,
        llh_dict) = load(f)
    mle_list.append(mle_dict)

fig, axes = plt.subplots(4, 2, figsize=(6, 9))
fig.tight_layout()
axes = axes.flatten()

for i, mle_dict in enumerate(mle_list):

    xVals=range(max(data_set_list[i])+1)
    no_vals = max(data_set_list[i]) + 1

    PoiLine=stats.poisson.pmf(xVals,mle_dict['poisson'])
    axes[i].plot(xVals,PoiLine,':s', label='Poisson')

    GeomLine=stats.geom.pmf(xVals,1/(mle_dict['geometric']+1),-1)
    axes[i].plot(xVals,GeomLine,'--v', label='Geometric')

    NegBinLine=stats.nbinom.pmf(xVals,
        mle_dict['negative binomial'][0]/mle_dict['negative binomial'][1],
        1/(mle_dict['negative binomial'][1]+1))
    axes[i].plot(xVals,NegBinLine,'-.x', label='Negative Binomial')

    ZIPLine=zip_pmf(xVals,mle_dict['zip'][0],mle_dict['zip'][1])
    axes[i].plot(xVals,ZIPLine,'^',linestyle=(0, (3, 5, 1, 5)), label='ZIP')

    if mle_dict['beta-Poisson'][2]>1e-4:
        BetaPoiLine=beta_poisson_pmf(xVals,
            mle_dict['beta-Poisson'][0],
            mle_dict['beta-Poisson'][1],
            1/mle_dict['beta-Poisson'][2])
        axes[i].plot(xVals,BetaPoiLine,'-o', label='Beta Poisson')

    counts,bins=np.histogram(data_set_list[i],no_vals)
    dist=counts/len(data_set_list[i])
    axes[i].bar(np.where(dist>0)[0],dist[dist>0], fill=False, label='Data')
    axes[i].axis([-0.5,no_vals+0.5,0,1.0])
    axes[i].set_aspect((no_vals+1))

    axes[i].set_xlabel('Secondary cases')
    axes[i].set_ylabel('Probability')

    axes[i].text(-4, 1, figlabels[i],
            fontsize=12,
            verticalalignment='top',
            fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

fig.savefig('plots/'+'fitted_distributions.png',format='png',bbox_inches='tight')

# x = np.linspace(1e-2,1, 100)
# y = stats.beta.pdf(x, np.mean(mle_list[i][i])*phi_mle, (1/N_inv_mle-np.mean(mle_list[i][i]))*phi_mle)
# fig,ax=plt.subplots(figsize=(5,5))
# plt.plot(x, y,'r-', lw=3, alpha=0.6, label='beta pdf')
# axes[i].axis([-0.01,1.01,0,1.01*np.max(y)])
# axes[i].set_aspect(1.02/(1.01*np.max(y)))
# plt.xlabel('Transmission probability')
# plt.ylabel('PDF')
# plt.xticks()
# plt.yticks()
# plt.title('Beta distribution underlying Gani and Leach data')
# fig.savefig('plots/'+'plague_beta_dist.png',format='png',bbox_inches='tight')
