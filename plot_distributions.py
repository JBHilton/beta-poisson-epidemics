'''This plots the output of the plague analysis.'''

import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)
from functions import beta_poisson_pmf, zip_pmf

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
for i, data_name in enumerate(data_name_list):
    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname,'rb') as f:
        (mle_dict,
            var_dict,
            superspread_dict,
            p0_dict,
            ci_dict,
            var_ci_dict,
            superspread_ci_dict,
            p0_ci_dict,
            llh_dict) = load(f)
    mle_list.append(mle_dict)
    ci_list.append(ci_dict)
    superspread_list.append(superspread_dict)
    superspread_ci_list.append(superspread_ci_dict)
    p0_list.append(p0_dict)
    p0_ci_list.append(p0_ci_dict)

fig, axes = plt.subplots(4, 2, figsize=(12, 18))
fig.tight_layout()
axes = axes.flatten()

for i, mle_dict in enumerate(mle_list):

    xMax = 6

    no_vals = max(data_set_list[i]) + 1
    counts,bins=np.histogram(data_set_list[i],no_vals)
    dist=counts/len(data_set_list[i])
    if max(data_set_list[i])>xMax:
        print(max(data_set_list[i]))
        print(xMax)
        dist_pos=np.where(dist[:xMax+1]>0)[0]
    else:
        dist_pos=np.where(dist>0)[0]
    # xVals = [str(i) for i in dist_pos]
    xVals = np.arange(xMax+1)
    axes[i].bar(dist_pos,dist[dist_pos], fill=False, label='Data')
    axes[i].axis([-0.5,xMax+0.5,0,1.0])
    axes[i].set_aspect((xMax+1))


    PoiLine=stats.poisson.pmf(xVals,mle_dict['poisson'])
    axes[i].plot(xVals,PoiLine,':s', label='Poisson')
    # axes[i].bar(xVals,PoiLine, label='Poisson')

    GeomLine=stats.geom.pmf(xVals,1/(mle_dict['geometric']+1),-1)
    axes[i].plot(xVals,GeomLine,'--v', label='Geometric')
    # axes[i].bar(xVals,GeomLine, label='Geometric')

    NegBinLine=stats.nbinom.pmf(xVals,
        mle_dict['negative binomial'][0]/mle_dict['negative binomial'][1],
        1/(mle_dict['negative binomial'][1]+1))
    axes[i].plot(xVals,NegBinLine,'-.x', label='Negative Binomial')
    # axes[i].bar(xVals,NegBinLine, label='Negative Binomial')

    ZIPLine=zip_pmf(xVals,mle_dict['zip'][0],mle_dict['zip'][1])
    axes[i].plot(xVals,ZIPLine,'^',linestyle=(0, (3, 5, 1, 5)), label='ZIP')
    # axes[i].bar(xVals,ZIPLine,linestyle=(0, (3, 5, 1, 5)), label='ZIP')

    if mle_dict['beta-Poisson'][2]>1e-4:
        BetaPoiLine=beta_poisson_pmf(xVals,
            mle_dict['beta-Poisson'][0],
            mle_dict['beta-Poisson'][1],
            1/mle_dict['beta-Poisson'][2])
        axes[i].plot(xVals,BetaPoiLine,'-o', label='Beta Poisson')
        # axes[i].bar(xVals,BetaPoiLine, label='Beta Poisson')

    axes[i].set_xlabel('Secondary cases')
    axes[i].set_ylabel('Probability')

    axes[i].text(-4, 1, figlabels[i],
            fontsize=12,
            verticalalignment='top',
            fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

axes[0].legend(loc='center left', bbox_to_anchor=(3., 0.5))

for fmt in formats:
    fig.savefig('plots/'+'fitted_distributions'+fmt,bbox_inches='tight')

fig, axes = plt.subplots(4, 2, figsize=(7.5, 15))
fig.tight_layout()
plt.subplots_adjust(hspace=0.3)
axes = axes.flatten()

for i, superspread_dict in enumerate(superspread_list):

    superspread_ci_dict = superspread_ci_list[i]

    label_list = ['Data',
                  'Beta-\nPoisson',
                  'Negative\nBinomial',
                  'Geometric',
                  'ZIP',
                  'Poisson']
    
    key_list = ['sample',
                'beta-Poisson',
                'negative binomial',
                'geometric',
                'zip',
                'poisson']
    
    bar_vals = np.array([
        superspread_dict[key] for key in key_list
    ])
    bar_errs = np.vstack((np.array([0,0]),
                         np.array([[superspread_dict[key] - superspread_ci_dict[key][0],  superspread_ci_dict[key][1] - superspread_dict[key]] for key in key_list[1:]]))
    )

    axes[i].bar(label_list, bar_vals, yerr=bar_errs.T)
    # axes[i].axis([-0.5, 4.5, 0, .1])
    axes[i].set_ylim([0, 0.18])
    axes[i].set_aspect(7/.18)
    axes[i].set_ylabel('Superspreading\n proportion')
    axes[i].set_xticklabels(label_list, rotation=45, ha='right')

    # axes[i].text(-4, 1, figlabels[i],
    #         fontsize=12,
    #         verticalalignment='top',
    #         fontfamily='serif',
    #         bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

for fmt in formats:
    fig.savefig('plots/'+'superspread_props'+fmt,bbox_inches='tight')

fig, axes = plt.subplots(4, 2, figsize=(7.5, 15))
fig.tight_layout()
plt.subplots_adjust(hspace=0.3)
axes = axes.flatten()

for i, p0_dict in enumerate(p0_list):

    p0_ci_dict = p0_ci_list[i]

    label_list = ['Data',
                  'Beta-\nPoisson',
                  'Negative\nBinomial',
                  'Geometric',
                  'ZIP',
                  'Poisson']
    
    key_list = ['sample',
                'beta-Poisson',
                'negative binomial',
                'geometric',
                'zip',
                'poisson']
    
    bar_vals = np.array([
        p0_dict[key] for key in key_list
    ])
    bar_errs = np.vstack((np.array([0,0]),
                         np.array([[p0_dict[key] - p0_ci_dict[key][0],  p0_ci_dict[key][1] - p0_dict[key]] for key in key_list[1:]]))
    )

    axes[i].bar(label_list, bar_vals, yerr=bar_errs.T)
    # axes[i].axis([-0.5, 4.5, 0, .1])
    axes[i].set_ylim([0, 1])
    axes[i].set_aspect(7/1)
    axes[i].set_ylabel('Superspreading\n proportion')
    axes[i].set_xticklabels(label_list, rotation=45, ha='right')

    # axes[i].text(-4, 1, figlabels[i],
    #         fontsize=12,
    #         verticalalignment='top',
    #         fontfamily='serif',
    #         bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

for fmt in formats:
    fig.savefig('plots/'+'p0_props'+fmt,bbox_inches='tight')

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
