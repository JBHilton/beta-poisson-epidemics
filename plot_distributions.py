''' In this script we plot the distributions in the schematic in Figure 1 of our manuscript.'''
#%%
from os import mkdir
from os.path import isdir
import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from datasets import (plague_data, mpox_data, nigeria_ebola_data,
    guinea_ebola_data, singapore_sars_data, sk_mers_data, sa_mers_data, noro_data)

x = np.linspace(1e-2,1, 100)
y = stats.beta.pdf(x, 1.325*0.58, (4-1.325)*0.58)
fig,ax=plt.subplots(figsize=(5,5))
plt.plot(x, y,'k-', lw=3, alpha=1., label='beta pdf')
ax.axis([-0.01,1.01,0,1.01*np.max(y)])
ax.set_aspect(1.02/(1.01*np.max(y)))
plt.xlabel('Transmission probability')
plt.ylabel('Probability density')
plt.xticks()
plt.yticks()
plt.title('Beta distributed infectivity')
fig.savefig('plots/underlying_beta.png',bbox_inches='tight')

x = np.arange(5)
y = stats.poisson.pmf(x, 1.325)
fig,ax=plt.subplots(figsize=(5,5))
plt.plot(x, y,'ko', ms=5, alpha=1., label='poisson pdf')
ax.axis([-0.5,4.5,0,1.1*np.max(y)])
ax.set_aspect(5./(1.1*np.max(y)))
plt.xlabel('Number of secondary cases')
plt.ylabel('Probability density')
plt.xticks()
plt.yticks()
plt.title('Poisson distributed secondary cases')
fig.savefig('plots/basic_underlying_poisson.png',bbox_inches='tight')

x = np.arange(11)
y = stats.poisson.pmf(x, 4.)
fig,ax=plt.subplots(figsize=(5,5))
plt.plot(x, y,'ko', ms=5, alpha=1., label='poisson pdf')
ax.axis([-0.5,10.5,0,1.1*np.max(y)])
ax.set_aspect(11./(1.1*np.max(y)))
plt.xlabel('Number of social contacts')
plt.ylabel('Probability density')
plt.xticks()
plt.yticks()
plt.title('Poisson distributed social contacts')
fig.savefig('plots/bp_underlying_poisson.png',bbox_inches='tight')

x = np.arange(6)
y = stats.poisson.pmf(x, 1.87)
fig,ax=plt.subplots(figsize=(5,5))
plt.plot(x, y,'ko', ms=5, alpha=1., label='poisson pdf')
ax.axis([-0.5,5.5,0,1.1*np.max(y)])
ax.set_aspect(6./(1.1*np.max(y)))
plt.xlabel('Number of secondary cases')
plt.ylabel('Probability density')
plt.xticks()
plt.yticks()
plt.title('Poisson distributed secondary cases')
fig.savefig('plots/zip_underlying_poisson.png',bbox_inches='tight')

x = np.linspace(1e-2,2.5, 100)
y = stats.gamma.pdf(x, 1.325/0.92, scale=0.92)
fig,ax=plt.subplots(figsize=(5,5))
plt.plot(x, y,'k-', lw=3, alpha=1., label='gamma pdf')
ax.axis([-0.05,2.55,0,1.1*np.max(y)])
ax.set_aspect(2.6/(1.1*np.max(y)))
plt.xlabel('Reproductive ratio')
plt.ylabel('Probability density')
plt.xticks()
plt.yticks()
plt.title('Gamma distributed individual reproductive ratio')
fig.savefig('plots/underlying_gamma.png',bbox_inches='tight')
# %%
