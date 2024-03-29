'''This plots the output of our main model fitting analysis.'''

from os import mkdir
from os.path import isdir
import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from datasets import (plague_data, mpox_data, nigeria_ebola_data,
    guinea_ebola_data, singapore_sars_data, sk_mers_data, sa_mers_data, noro_data)
from functions import beta_poisson_pmf, zip_pmf

if isdir('outputs/mles') is True:
    fname_root = 'outputs/mles/'
else:
    fname_root = 'reference-outputs/mles/'

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
for i, data_key_list in enumerate(data_key_list):
    fname = fname_root+data_key_list+'_results.pkl'
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

label_list = ['$\\hat{\\lambda}$',
                '$\\hat{\\Phi}$',
                '$\\hat{\\nu}$']

fig, axes = plt.subplots(1, 3, figsize=(20, 5))
fig.tight_layout()
axes = axes.flatten()
    
for idx in range(3):
    
    bar_vals = np.array([
        mle_dict['beta-Poisson'][idx] for mle_dict in mle_list
    ])
    bar_errs = np.array([[mle_list[i]['beta-Poisson'][idx] - ci_list[i]['beta-Poisson'][idx][0], 
                          ci_list[i]['beta-Poisson'][idx][1] - mle_list[i]['beta-Poisson'][idx]] for i in range(len(mle_list))])

    axes[idx].errorbar(data_name_list,
                       bar_vals,
                       yerr=bar_errs.T,
                       lw=0,
                       elinewidth=2,
                       marker='o',
                       ms=8,
                       color='steelblue',
                       capsize=8)
    axes[idx].set_ylabel(label_list[idx], rotation=0, labelpad=20)
    axes[idx].set_xticklabels(data_name_list, rotation=45, ha='right')


axes[1].set_yscale('log')
axes[2].set_yscale('log')
axes[0].set_ylim([0, 2.5])
axes[0].set_aspect(8 / 2.5)
axes[1].set_ylim([1e-7, 20])
axes[1].set_aspect(8 / (np.log10(20)+7))
axes[2].set_ylim([1e-5, 20])
axes[2].set_aspect(8 / (np.log10(20)+5))

axes[0].text(-3, 2.5, 'a)',
        fontsize=16,
        verticalalignment='top',
        fontfamily='serif',
        bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
axes[1].text(-3, 20, 'b)',
        fontsize=16,
        verticalalignment='top',
        fontfamily='serif',
        bbox=dict(facecolor='1', edgecolor='none', pad=3.0))
axes[2].text(-3, 20, 'c)',
        fontsize=16,
        verticalalignment='top',
        fontfamily='serif',
        bbox=dict(facecolor='1', edgecolor='none', pad=3.0))


for fmt in formats:
    fig.savefig('plots/'+'bp_mles'+fmt,bbox_inches='tight')

fig, axes = plt.subplots(4, 2, figsize=(15, 35))
fig.tight_layout()
plt.subplots_adjust(wspace=0.3)
axes = axes.flatten()

for i, mle_dict in enumerate(mle_list):

    xMax = 6

    no_vals = max(data_set_list[i]) + 1
    counts,bins=np.histogram(data_set_list[i],no_vals)
    dist=counts/len(data_set_list[i])
    if max(data_set_list[i])>xMax:
        dist_pos=np.where(dist[:xMax+1]>0)[0]
    else:
        dist_pos=np.where(dist>0)[0]
    # xVals = [str(i) for i in dist_pos]
    xVals = np.arange(xMax+1)
    axes[i].bar(dist_pos,dist[dist_pos], fill=False, label='Data')
    axes[i].axis([-0.5,xMax+0.5,0,1.0])
    axes[i].set_aspect((xMax+1))


    PoiLine=stats.poisson.pmf(xVals,mle_dict['poisson'])
    axes[i].plot(xVals,PoiLine,':s', label='Poisson', lw=2, ms=12)
    # axes[i].bar(xVals,PoiLine, label='Poisson')

    GeomLine=stats.geom.pmf(xVals,1/(mle_dict['geometric']+1),-1)
    axes[i].plot(xVals,GeomLine,'--v', label='Geometric', lw=2, ms=12)
    # axes[i].bar(xVals,GeomLine, label='Geometric')

    NegBinLine=stats.nbinom.pmf(xVals,
        mle_dict['negative binomial'][0]/mle_dict['negative binomial'][1],
        1/(mle_dict['negative binomial'][1]+1))
    axes[i].plot(xVals,NegBinLine,'-.x', label='Negative Binomial', lw=2, ms=12)
    # axes[i].bar(xVals,NegBinLine, label='Negative Binomial')

    ZIPLine=zip_pmf(xVals,mle_dict['zip'][0],mle_dict['zip'][1])
    axes[i].plot(xVals,ZIPLine,'^',linestyle=(0, (3, 5, 1, 5)), label='ZIP', lw=2, ms=12)
    # axes[i].bar(xVals,ZIPLine,linestyle=(0, (3, 5, 1, 5)), label='ZIP')

    if mle_dict['beta-Poisson'][2]>1e-4:
        BetaPoiLine=beta_poisson_pmf(xVals,
            mle_dict['beta-Poisson'][0],
            mle_dict['beta-Poisson'][1],
            1/mle_dict['beta-Poisson'][2])
        axes[i].plot(xVals,BetaPoiLine,'-o', label='Beta Poisson', lw=2, ms=12)
        # axes[i].bar(xVals,BetaPoiLine, label='Beta Poisson')

    axes[i].set_xlabel('Secondary cases')
    axes[i].set_ylabel('Probability')
    axes[i].set_title(data_name_list[i])

    axes[i].text(-2, 1, figlabels[i],
            verticalalignment='top',
            fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

axes[0].legend(loc='center left', bbox_to_anchor=(2.4, 0.5))

for fmt in formats:
    fig.savefig('plots/'+'fitted_distributions'+fmt,bbox_inches='tight')

fig, axes = plt.subplots(4, 2, figsize=(15, 30))
fig.tight_layout()
plt.subplots_adjust(hspace=.3)
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

    axes[i].errorbar(label_list,
                       bar_vals,
                       yerr=bar_errs.T,
                       lw=0,
                       elinewidth=2,
                       marker='_',
                       ms=12,
                       color='steelblue',
                       capsize=8)
    # axes[i].axis([-0.5, 4.5, 0, .1])
    axes[i].set_ylim([0, 0.18])
    axes[i].set_aspect(7/.18)
    axes[i].set_ylabel('Superspreading\n proportion')
    axes[i].set_xticklabels(label_list, rotation=45, ha='right')
    axes[i].set_title(data_name_list[i])

    axes[i].text(-2, .18, figlabels[i],
            verticalalignment='top',
            fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

for fmt in formats:
    fig.savefig('plots/'+'superspread_props'+fmt,bbox_inches='tight')

fig, axes = plt.subplots(4, 2, figsize=(15, 30))
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

    axes[i].errorbar(label_list,
                       bar_vals,
                       yerr=bar_errs.T,
                       lw=0,
                       elinewidth=2,
                       marker='_',
                       ms=12,
                       color='steelblue',
                       capsize=8)
    # axes[i].axis([-0.5, 4.5, 0, .1])
    axes[i].set_ylim([0, 1])
    axes[i].set_aspect(7/1)
    axes[i].set_ylabel('P[0]')
    axes[i].set_xticklabels(label_list, rotation=45, ha='right')
    axes[i].set_title(data_name_list[i])

    axes[i].text(-2, 1, figlabels[i],
            verticalalignment='top',
            fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

for fmt in formats:
    fig.savefig('plots/'+'p0_bars'+fmt,bbox_inches='tight')

fig, axes = plt.subplots(4, 2, figsize=(15, 30))
fig.tight_layout()
plt.subplots_adjust(hspace=.3)
axes = axes.flatten()

for i, od_dict in enumerate(od_list):

    od_ci_dict = od_ci_list[i]

    label_list = ['Data',
                  'Beta-\nPoisson',
                  'Negative\nBinomial',
                  'Geometric',
                  'ZIP']
    
    key_list = ['sample',
                'beta-Poisson',
                'negative binomial',
                'geometric',
                'zip']
    
    bar_vals = np.array([
        od_dict[key] for key in key_list
    ])
    bar_errs = np.vstack((np.array([0,0]),
                         np.array([[od_dict[key] - od_ci_dict[key][0],  od_ci_dict[key][1] - od_dict[key]] for key in key_list[1:]]))
    )

    y_max = 5 * np.ceil((bar_vals + bar_errs[:, 1]).max()/5)

    axes[i].errorbar(label_list,
                       bar_vals,
                       yerr=bar_errs.T,
                       lw=0,
                       elinewidth=2,
                       marker='_',
                       ms=12,
                       color='steelblue',
                       capsize=8)
    axes[i].set_ylim([0, y_max])
    axes[i].set_aspect(6/y_max)
    axes[i].set_ylabel('Overdispersion')
    axes[i].set_xticklabels(label_list, rotation=45, ha='right')
    axes[i].set_title(data_name_list[i])

    axes[i].text(-2, y_max, figlabels[i],
            verticalalignment='top',
            fontfamily='serif',
            bbox=dict(facecolor='1', edgecolor='none', pad=3.0))

for fmt in formats:
    fig.savefig('plots/'+'od_bars'+fmt,bbox_inches='tight')