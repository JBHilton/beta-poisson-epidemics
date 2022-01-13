'''This plots the output of the plague analysis.'''

import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from datasets import plague_data
from functions import beta_poisson_pmf, zip_pmf


with open('outputs/plague/results.pkl','rb') as f:
    (theta_mle,
    lmbd_mle,
    phi_mle,
    N_inv_mle,
    sigma_mle) = load(f)

fig, ax=plt.subplots(figsize=(5,5))

xVals=range(max(plague_data)+1)

PoiLine=stats.poisson.pmf(xVals,np.mean(plague_data))
ax.plot(xVals,PoiLine,':s', label='Poisson')
GeomLine=stats.geom.pmf(xVals,1/(np.mean(plague_data)+1),-1)
ax.plot(xVals,GeomLine,'--v', label='Geometric')
NegBinLine=stats.nbinom.pmf(xVals,np.mean(plague_data)/theta_mle,1/(theta_mle+1))
ax.plot(xVals,NegBinLine,'-.x', label='Negative Binomial')
ZIPLine=zip_pmf(xVals,lmbd_mle,sigma_mle)
ax.plot(xVals,ZIPLine,'^',linestyle=(0, (3, 5, 1, 5)), label='ZIP')
BetaPoiLine=beta_poisson_pmf(xVals,np.mean(plague_data),phi_mle,1/N_inv_mle)
ax.plot(xVals,BetaPoiLine,'-o', label='Beta Poisson')
counts,bins=np.histogram(plague_data,7)
dist=counts/len(plague_data)
ax.bar(np.where(dist>0)[0],dist[dist>0], fill=False, label='Data')
ax.axis([-0.5,6.5,0,0.5])
ax.set_aspect(7/0.5)

ax.legend()
plt.xlabel('Secondary cases')
plt.ylabel('Probability')
plt.xticks()
plt.yticks()
plt.title('Plague (Gani and Leach 2004)')
fig.savefig('plague_fit.png',format='png',bbox_inches='tight')

x = np.linspace(1e-2,1, 100)
y = stats.beta.pdf(x, np.mean(plague_data)*phi_mle, (1/N_inv_mle-np.mean(plague_data))*phi_mle)
fig,ax=plt.subplots(figsize=(5,5))
plt.plot(x, y,'r-', lw=3, alpha=0.6, label='beta pdf')
ax.axis([-0.01,1.01,0,1.01*np.max(y)])
ax.set_aspect(1.02/(1.01*np.max(y)))
plt.xlabel('Transmission probability')
plt.ylabel('PDF')
plt.xticks()
plt.yticks()
plt.title('Beta distribution underlying Gani and Leach data')
fig.savefig('plague_beta_dist.png',format='png',bbox_inches='tight')
