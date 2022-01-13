'''This script contains all of the functions used in our study'''

from __future__ import print_function
import math
import numpy as np
import scipy as sp
from scipy import special as spsp
from scipy import sparse as sparse
from scipy import stats as stats
from scipy import optimize as opt
from scipy import interpolate as interpolate
import time
import mpmath
import random

def poisson_loglh(data,lmbd):
    '''
    Calculate log likelihood of Poisson parameter lambda given data.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd : float
            estimated Poisson parameter

    Returns
    -------
        llh : float
            log likelihood of lmbd given data
    '''
    llh=0
    for x in data:
        llh+=np.log(stats.poisson.pmf(x,lmbd))
    return llh

def geo_loglh(data,lmbd):
    '''
    Calculate log likelihood of geometric parameter lambda given data.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd : float
            estimated geometric parameter

    Returns
    -------
        llh : float
            log likelihood of lmbd given data
    '''
    llh=0
    for x in data:
        llh+=np.log(stats.geom.pmf(x,1/(lmbd+1),-1))
    return llh

def neg_bin_loglh_theta(data,lmbd,theta):
    '''
    Calculate log likelihood of negative binomial parameters given data.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd : float
            estimated mean of negative binomial distribution
        theta : float
            estimated overdispersion of negative binomial distribution

    Returns
    -------
        llh : float
            log likelihood of lmbd and theta given data
    '''
    llh=0
    for x in data:
        llh+=np.log(stats.nbinom.pmf(x,lmbd/theta,1/(theta+1)))
    return llh

def get_theta_mle(data,theta_0):
    '''
    Calculate maximum likelihood estimate of negative binomial overdispersion
    parameter theta given sample data.

    Parameters
    ----------
        data : list
            sample dataset
        theta_0 : float
            initial estimate of overdispersion parameter

    Returns
    -------
        : float
            maximum likelihood estimate of overdispersion parameter
    '''
    def f(theta):
        lmbd=np.mean(data)
        return -neg_bin_loglh_theta(data,lmbd,theta)
    mle=sp.optimize.minimize(f,[theta_0],bounds=((1e-6,50),))
    return mle.x[0]

def beta_poisson_pmf(x,lmbd,Phi,N):
    '''
    Evaluate the probability mass function for beta-Poisson distribution.

    Parameters
    ----------
        x : int or array
            point(s) at which to evaluate function
        lmbd : float
        Phi : float
        N : float

    Returns
    -------
        P : float or array
            probability of each point in x

    '''
    if type(x)==int:
        P=spsp.hyp1f1(x+Phi*lmbd,x+Phi*N,-N)
        for n in range(1,x+1): # This loop gives us the N^x/gamma(x+1 term)
            P=(N/n)*P
        for m in range(x): # This loop gives us the term with the two gamma functions in numerator and denominator
            P=((m+Phi*lmbd)/(m+Phi*N))*P
    else:
        P=[]
        for i in range(0,len(x)):
            p=spsp.hyp1f1(x[i]+Phi*lmbd,x[i]+Phi*N,-N)
            for n in range(1,x[i]+1): # This loop gives us the N^x/gamma(x+1 term)
                p=(N/n)*p
            for m in range(x[i]): # This loop gives us the term with the two gamma functions in numerator and denominator
                p=((m+Phi*lmbd)/(m+Phi*N))*p
            P=P+[p]
    return P

hyp1f1_alt=np.frompyfunc(mpmath.hyp1f1,3,1)
def beta_poisson_loglh(data,lmbd,phi,N):
    '''
    Calculate log likelihood of beta-Poisson parameters given data.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd : float
        phi : float
        N : float

    Returns
    -------
        llh : float
            log likelihood of parameters given data
    '''
    llh=0
    for x in data:
        llh+=x*np.log(N)-np.real(spsp.loggamma(x+1))+np.real(spsp.loggamma(phi*N))+np.real(spsp.loggamma(x+phi*lmbd))-np.real(spsp.loggamma(x+phi*N))-np.real(spsp.loggamma(phi*lmbd))
        if x+phi*N<50:
            llh+=np.log(spsp.hyp1f1(x+phi*lmbd,x+phi*N,-N))
        else:
            llh+=np.log(float(hyp1f1_alt(x+phi*lmbd,x+phi*N,-N)))
    return llh

def neg_bin_loglh(data,lmbd,phi):
    '''
    Calculate log likelihood of negative binomial parameters given data, with
    negative binomial parameterised with phi rather than theta.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd : float
            estimated mean of negative binomial distribution
        phi : float

    Returns
    -------
        llh : float
            log likelihood of lmbd and theta given data
    '''
    llh=0
    for x in data:
        llh+=np.log(stats.nbinom.pmf(x,lmbd*phi,phi/(phi+1)))
    return llh

def get_phi_and_N_mles(data,phi_0,N_0):
    '''
    Calculate maximum likelihood estimates of beta-Poisson parameters Phi and N.

    Parameters
    ----------
        data : list
            sample dataset
        phi_0 : float
            initial estimate of Phi
        N_0 : float
            initial estimate of N

    Returns
    -------
        : float
            maximum likelihood estimate of Phi
        : float
            maximum likelihood estimate of N
    '''
    def f(params):
        lmbd=np.mean(data)
        phi=params[0]
        if params[1]>0.1e-3:
            N=1/params[1]
            return -beta_poisson_loglh(data,lmbd,phi,N)
        else:
            return -neg_bin_loglh(data,lmbd,phi)
    mle=sp.optimize.minimize(f,[phi_0,N_0],bounds=((1e-6,50),(0,1/np.mean(data))))
    if mle.x[1]<0:
        mle.x[1]=0
    return mle.x[0],mle.x[1]

def zip_pmf(x,lmbd,sigma):
    '''
    Evaluate the probability mass function for zero-inflated Poisson
    distribution.

    Parameters
    ----------
        x : int or array
            point(s) at which to evaluate function
        lmbd : float
            mean of Poisson component
        sigma : float
            degree of zero inflation

    Returns
    -------
        P : float or array
            probability of each point in x

    '''
    if type(x)==int:
        return sigma*(x==0)+(1-sigma)*stats.poisson.pmf(x,lmbd)
    else:
        return sigma*np.equal(x,np.zeros(len(x)))+(1-sigma)*stats.poisson.pmf(x,lmbd)

def zip_loglh(data,lmbd,sigma):
    '''
    Calculate log likelihood of zero-inflated Poisson parameters given data.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd : float
            mean of Poisson component
        sigma : float
            degree of zero inflation
    Returns
    -------
        llh : float
            log likelihood of lmbd and sigma given data
    '''
    llh=0
    for x in data:
        if x==0:
            llh+=np.log(sigma+(1-sigma)*np.exp(-lmbd))
        else:
            llh+=np.log(1-sigma)+np.log(stats.poisson.pmf(x,lmbd))
    return llh

def get_zip_mles(data,lmbd_0,sigma_0):
    '''
    Calculate maximum likelihood estimates of ZIP parameters lambda and sigma.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd_0 : float
            initial estimate of lambda
        sigma_0 : float
            initial estimate of sigma

    Returns
    -------
        : float
            maximum likelihood estimate of lambda
        : float
            maximum likelihood estimate of sigma
    '''
    def f(params):
        lmbd=params[0]
        sigma=params[1]
        return -zip_loglh(data,lmbd,sigma)
    mle=sp.optimize.minimize(f,[lmbd_0,sigma_0],bounds=((np.mean(data),50),(0,1-1e-6)))
    return mle.x[0],mle.x[1]

def beta_poisson_pgf(s,lmbd,phi,N):
    '''
    Probability generating function of the beta-Poisson distribution.

    Parameters
    ----------
        s : float
            point at which to evaluate PGF
        lmbd : float
        Phi : float
        N : float

    Returns
    -------
        G : float or array
            PGF evaluated at s
    '''
    if np.size(np.asarray(s))==1:
        G=spsp.hyp1f1(lmbd*phi,N*phi,N*(s-1));
    else:
        G=[]
        for i in range(0,np.size(np.asarray(s))):
            G=G+[spsp.hyp1f1(lmbd*phi,N*phi,N*(s[i]-1))]
    return G

def poisson_pgf(s,lmbd):
    '''
    Probability generating function of the Poisson distribution.

    Parameters
    ----------
        s : float
            point at which to evaluate PGF
        lmbd : float

    Returns
    -------
        G : float or array
            PGF evaluated at s
    '''
    if np.size(np.asarray(s))==1:
        G=np.exp(lmbd*(s-1))
    else:
        G=[]
        for i in range(0,np.size(np.asarray(s))):
            G=G+[np.exp(lmbd*(s[i]-1))]
    return G

def geom_pgf(s,lmbd):
    '''
    Probability generating function of the geometric distribution.

    Parameters
    ----------
        s : float
            point at which to evaluate PGF
        lmbd : float

    Returns
    -------
        G : float or array
            PGF evaluated at s
    '''
    if np.size(np.asarray(s))==1:
        G=1/(lmbd+1-lmbd*s)
    else:
        G=[]
        for i in range(0,np.size(np.asarray(s))):
            G=G+[1/(lmbd+1-lmbd*s[i])]
    return G

def neg_bin_pgf(s,lmbd,theta):
    '''
    Probability generating function of the negative binomial distribution.

    Parameters
    ----------
        s : float
            point at which to evaluate PGF
        lmbd : float
        theta : float

    Returns
    -------
        G : float or array
            PGF evaluated at s
    '''
    if np.size(np.asarray(s))==1:
        G=(theta+1-s*theta)**(-lmbd/theta)
    else:
        G=[]
        for i in range(0,np.size(np.asarray(s))):
            G=G+[(theta+1-s[i]*theta)**(-lmbd/theta)]
    return G

def zip_pgf(s,lmbd,sigma):
    '''
    Probability generating function of the zero-inflated Poisson distribution.

    Parameters
    ----------
        s : float
            point at which to evaluate PGF
        lmbd : float
        sigma : float

    Returns
    -------
        G : float or array
            PGF evaluated at s
    '''
    if np.size(np.asarray(s))==1:
        G=sigma+(1-sigma)*np.exp(lmbd*(s-1))
    else:
        G=[]
        for i in range(0,np.size(np.asarray(s))):
            G=G+[sigma+(1-sigma)*np.exp(lmbd*(s[i]-1))]
    return G

def beta_poisson_extinction_prob( lmbd,phi,N ):
    '''
    Calculate the probability that the beta-Poisson branching process becomes
    extinct.

    Parameters
    ----------
        lmbd : float
        phi : float
        N : float

    Returns
    -------
        q : float
            extinction probability
    '''
    if lmbd<=1:
        return 1
    else:
        def f(s):
            return beta_poisson_pgf(s,lmbd,phi,N)-s
        q=opt.brentq(
            f, 0, 1-1e-4);
        return q

def poisson_extinction_prob( lmbd ):
    '''
    Calculate the probability that the Poisson branching process becomes
    extinct.

    Parameters
    ----------
        lmbd : float

    Returns
    -------
        q : float
            extinction probability
    '''
    if lmbd<=1:
        return 1
    else:
        def f(s):
            return poisson_pgf(s,lmbd)-s
    q=opt.brentq(
        f, 0, 1-1e-6)
    return q

def geom_extinction_prob( lmbd ):
    '''
    Calculate the probability that the geometric branching process becomes
    extinct.

    Parameters
    ----------
        lmbd : float

    Returns
    -------
        q : float
            extinction probability
    '''
    if lmbd<=1:
        return 1
    else:
        def f(s):
            return geom_pgf(s,lmbd)-s
    q=opt.brentq(
        f, 0, 1-1e-6)
    return q

def neg_bin_extinction_prob( lmbd,theta ):
    '''
    Calculate the probability that the negative binomial branching process
    becomes extinct.

    Parameters
    ----------
        lmbd : float
        theta : float

    Returns
    -------
        q : float
            extinction probability
    '''
    if lmbd<=1:
        return 1
    else:
        def f(s):
            return neg_bin_pgf(s,lmbd,theta)-s
        q=opt.brentq(
            f, 0, 1-1e-6)
        return q

def zip_extinction_prob( lmbd,sigma ):
    '''
    Calculate the probability that the zero-inflated Poisson branching process
    becomes extinct.

    Parameters
    ----------
        lmbd : float
        sigma : float

    Returns
    -------
        q : float
            extinction probability
    '''
    if (1-sigma)*lmbd<=1:
        return 1
    else:
        def f(s):
            return zip_pgf(s,lmbd,sigma)-s
        q=opt.brentq(
            f, 0, 1-1e-6)
        return q

def empirical_loglh(data):
    '''
    Calculate upper bound on log likelihood by using empirical distribution
    based on observed frequencies in data.

    Parameters
    ----------
        data : list
            sample data

    Returns
    -------
        llh : float
            log likelihood
    '''
    counts,bins=np.histogram(data,max(data)+1)
    dist=counts/len(data)
    llh=0
    for x in data:
        llh+=np.log(dist[x])
    return llh

def get_lambda_and_phi_mles(data,lmbd_0,phi_0,N_emp):
    '''
    Calculate maximum likelihood estimates of beta-Poisson parameters lambda and
    Phi given empirical estimate of contact parameter N.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd_0 : float
            initial estimate of lambda
        phi_0 : float
            initial estimate of Phi
        N_emp : float
            empirical estimate of N

    Returns
    -------
        : float
            maximum likelihood estimate of lambda
        : float
            maximum likelihood estimate of Phi
    '''
    def f(params):
        lmbd=params[0]
        phi=params[1]
        return -beta_poisson_loglh(data,lmbd,phi,N_emp)
    mle=sp.optimize.minimize(f,[lmbd_0,phi_0],bounds=((1e-6,10),(1e-6,50)))
    return mle.x[0],mle.x[1]

def generate_mle_dict(data,
                      theta_0,
                      phi_0,
                      N_0,
                      lmbd_0,
                      sigma_0):
    '''
    Calculate maximum likelihood estimates of parameters for each offspring
    distribution and output them as a dictionary.

    Parameters
    ----------
        data : list
            sample data to fit to
        theta_0 : float
            initial estimate of negative binomial overdispersion parameter
        phi_0 : float
            initial estimate of beta-Poisson parameter Phi
        N_0 : float
            initial estimate of beta-Poisson parameter Phi
        lmbd_0 : float
            initial estimate of ZIP baseline lambda
        sigma_0 : float
            initial estimate of ZIP inflation parameter sigma

    Returns
    -------
        mle_dict : dictionary
            dictionary containing maximum likelihood estimates of parameters for
            each model.
    '''
    theta_mle=get_theta_mle(data,1.5)
    phi_mle,N_inv_mle=get_phi_and_N_mles(data,0.5,1/2)
    lmbd_mle,sigma_mle=get_zip_mles(data,1.5,0.5)

    mle_dict = {
        'poisson' : np.mean(data),
        'geometric' : np.mean(data),
        'negative binomial' : [np.mean(data), theta_mle],
        'zip' : [lmbd_mle, sigma_mle],
        'beta-Poisson' : [np.mean(data), phi_mle, N_inv_mle]
    }

    return mle_dict

def poisson_mle_grid(data, interval, points):
    '''
    Calculate a confidence interval for the MLE of the Poisson distribution
    given some data using a grid calculation.

    Parameters
    ----------
        data : list
            sample dataset
        interval : list
            list containing the lower and upper bounds on which to perform the
            grid calculation
        points : int
            number of points to use in the grid calculation
    Returns
    -------
        lmbd : array
            value of lambda over grid
        llh : array
            log likelihood of lambda values over grid
        mle : float
            maximum likelihood estimate of lambda
        ci : list
            95% confidence interval for lambda given data
    '''
    lmbd=np.linspace(interval[0],interval[1],points)
    llh=np.zeros(points)
    for x in data:
        llh=llh+np.log(stats.poisson.pmf(x,lmbd))
    mle_loc=np.argmax(llh)
    mle=lmbd[mle_loc]
    lh_normed=np.exp(llh)/np.sum(np.exp(llh))
    current_max=lh_normed[mle_loc]
    interval_weight=current_max
    while interval_weight<0.95:
        max_loc=np.argmax(lh_normed[np.where(lh_normed<current_max)])
        current_max=lh_normed[np.where(lh_normed<current_max)][max_loc]
        interval_weight+=current_max
    ci=[np.min(lmbd[np.where(lh_normed>=current_max)[0]]),np.max(lmbd[np.where(lh_normed>=current_max)[0]])]
    return lmbd,llh,mle,ci

def geometric_mle_grid(data,interval,points):
    '''
    Calculate a confidence interval for the MLE of the geometric distribution
    given some data using a grid calculation.

    Parameters
    ----------
        data : list
            sample dataset
        interval : list
            list containing the lower and upper bounds on which to perform the
            grid calculation
        points : int
            number of points to use in the grid calculation
    Returns
    -------
        lmbd : array
            value of lambda over grid
        llh : array
            log likelihood of lambda values over grid
        mle : float
            maximum likelihood estimate of lambda
        ci : list
            95% confidence interval for lambda given data
    '''
    lmbd=np.linspace(interval[0],interval[1],points)
    llh=np.zeros(points)
    for x in data:
        llh=llh+np.log(stats.geom.pmf(x,1/(lmbd+1),-1))
    mle_loc=np.argmax(llh)
    mle=lmbd[mle_loc]
    lh_normed=np.exp(llh)/np.sum(np.exp(llh))
    current_max=lh_normed[mle_loc]
    interval_weight=current_max
    while interval_weight<0.95:
        max_loc=np.argmax(lh_normed[np.where(lh_normed<current_max)])
        current_max=lh_normed[np.where(lh_normed<current_max)][max_loc]
        interval_weight+=current_max
    ci=[np.min(lmbd[np.where(lh_normed>=current_max)[0]]),np.max(lmbd[np.where(lh_normed>=current_max)[0]])]
    return lmbd,llh,mle,ci

def zip_mle_grid(data,lmbd_interval,sigma_interval,lmbd_points,sigma_points):
    '''
    Calculate confidence intervals for the MLEs of the ZIP distribution
    parameters given some data using a grid calculation.

    Parameters
    ----------
        data : list
            sample dataset
        lmbd_interval : list
            list containing the lower and upper bounds on lambda
        sigma_interval : list
            list containing the lower and upper bounds on sigma
        lmbd_points : int
            number of grid points in the lambda direction
        sigma_points : int
            number of grid points in the sigma direction
    Returns
    -------
        lmbd_grid : array
            value of lambda over grid
        sigma_grid : array
            value of sigma over grid
        llh : array
            log likelihood of parameter values over grid
        lmbd_mle : float
            maximum likelihood estimate of lambda
        lmbd_ci : list
            95% confidence interval for lambda given data
        sigma_mle : float
            maximum likelihood estimate of sigma
        sigma_ci : list
            95% confidence interval for sigma given data
    '''
    lmbd=np.linspace(lmbd_interval[0],lmbd_interval[1],lmbd_points)
    sigma=np.linspace(sigma_interval[0],sigma_interval[1],sigma_points)

    lmbdgrid,sigmagrid=np.meshgrid(lmbd,sigma)
    llh=np.zeros(lmbdgrid.shape)
    for x in data:
        llh=llh+zip_log_lh(x,lmbdgrid,sigmagrid)
    mle_loc=np.unravel_index(np.argmax(llh),llh.shape)
    lmbdmle=lmbdgrid[mle_loc]
    sigmamle=sigmagrid[mle_loc]
    lh_normed=np.exp(llh)/np.sum(np.exp(llh))
    current_max=lh_normed[mle_loc]
    interval_weight=current_max
    while interval_weight<0.95:
        max_loc=np.unravel_index(np.argmax(lh_normed),lh_normed.shape)
        current_max=lh_normed[max_loc]
        lh_normed[max_loc]=0
        interval_weight+=current_max
    lh_normed=np.exp(llh)/np.sum(np.exp(llh))
    lmbdci=[np.min(lmbdgrid[np.where(lh_normed>=current_max)]),np.max(lmbdgrid[np.where(lh_normed>=current_max)])]
    sigmaci=[np.min(sigmagrid[np.where(lh_normed>=current_max)]),np.max(sigmagrid[np.where(lh_normed>=current_max)])]
    return lmbdgrid,sigmagrid,llh,lmbdmle,sigmamle,lmbdci,sigmaci

def neg_bin_bootstrap(data,no_samples,theta_0):
    '''
    Calculate confidence intervals for negative binomial parameter MLEs using
    bootstrapping.

    Parameters
    ----------
        data : list
            sample dataset
        no_samples : int
            number of samples to use in bootstrap
        theta_0 : float
            starting value of theta to use in MLE calculations

    Returns
    -------
        lmbd_mle : float
            maximum likelihood estimate of lambda
        lmbd_ci : list
            list containing lower and upper 95% confidence intervals for MLE of
            lambda
        lmbd_samples : array
            complete set of lambda estimates generated during sampling process
        theta_mle : float
            maximum likelihood estimate of theta
        theta_ci : list
            list containing lower and upper 95% confidence intervals for MLE of
            theta
        theta_samples : array
            complete set of theta estimates generated during sampling process
        var_mle : float
            maximum likelihood estimate of variance
        var_ci : list
            list containing lower and upper 95% confidence intervals for MLE of
            variance
        var_samples : array
            complete set of var estimates generated during sampling process
            calculated from lambda and theta samples
    '''

    sample_size=np.size(data)

    lmbd_mle=np.mean(data)
    theta_mle=get_theta_mle(data,theta_0)
    var_mle=lmbd_mle*(1+theta_mle)
    lmbd_samples=np.zeros(no_samples)
    theta_samples=np.zeros(no_samples)
    print('Now calculating',no_samples,'bootstrap samples.')
    start_time=time.time()

    for i in range(no_samples):
        data_now=random.choices(data,k=sample_size)
        lmbd_samples[i]=np.mean(data_now)
        theta_samples[i]=get_theta_mle(data_now,theta_0)
        if ((i+1)%100)==0:
            print('Sample',i+1,'of',no_samples,'completed.',time.time()-start_time,'seconds elapsed, approximately',(no_samples-i-1)*(time.time()-start_time)/(i+1),'remaining.')

    var_samples=np.multiply(lmbd_samples,(1+theta_samples))

    lmbd_ci=[np.percentile(lmbd_samples,2.5),np.percentile(lmbd_samples,97.5)]
    theta_ci=[np.percentile(theta_samples,2.5),np.percentile(theta_samples,97.5)]
    var_ci=[np.percentile(var_samples,2.5),np.percentile(var_samples,97.5)]

    return lmbd_mle,lmbd_ci,lmbd_samples,theta_mle,theta_ci,theta_samples,var_mle,var_ci,var_samples

def beta_poisson_bootstrap(data,no_samples,phi_0,N_0):
    '''
    Calculate confidence intervals for beta-Poisson parameter MLEs using
    bootstrapping.

    Parameters
    ----------
        data : list
            sample dataset
        no_samples : int
            number of samples to use in bootstrap
        phi_0 : float
            starting value of Phi to use in MLE calculations
        N_0 : float
            starting value of N to use in MLE calculations

    Returns
    -------
        lmbd_mle : float
            maximum likelihood estimate of lambda
        lmbd_ci : list
            list containing lower and upper 95% confidence intervals for MLE of
            lambda
        lmbd_samples : array
            complete set of lambda estimates generated during sampling process
        phi_mle : float
            maximum likelihood estimate of phi
        phi_ci : list
            list containing lower and upper 95% confidence intervals for MLE of
            phi
        phi_samples : array
            complete set of phi estimates generated during sampling process
        N_inv_mle : float
            maximum likelihood estimate of 1/N
        N_inv_ci : list
            list containing lower and upper 95% confidence intervals for MLE of
            1/N
        N_inv_samples : array
            complete set of N_inv estimates generated during sampling process
            calculated from lambda and theta samples
    '''

    sample_size=np.size(data)

    lmbd_mle=np.mean(data)
    phi_mle,N_inv_mle=get_phi_and_N_mles(data,phi_0,N_0)
    var_mle=lmbd_mle*(1+(1-lmbd_mle*N_inv_mle)/(phi_mle+N_inv_mle))
    lmbd_samples=np.zeros(no_samples)
    phi_samples=np.zeros(no_samples)
    N_inv_samples=np.zeros(no_samples)
    print('Now calculating',no_samples,'bootstrap samples.')
    start_time=time.time()

    for i in range(no_samples):
        data_now=random.choices(data,k=sample_size)
        lmbd_samples[i]=np.mean(data_now)
        phi_samples[i],N_inv_samples[i]=get_phi_and_N_mles(data_now,phi_0,N_0)
        if ((i+1)%100)==0:
            print('Sample',i+1,'of',no_samples,'completed.',time.time()-start_time,'seconds elapsed, approximately',(no_samples-i-1)*(time.time()-start_time)/(i+1),'remaining.')

    sample_array=np.zeros((no_samples,3))
    sample_array[:,1]=lmbd_samples # np.histogram does things in a slightly unintuitive order - this makes it come out right
    sample_array[:,0]=phi_samples
    sample_array[:,2]=N_inv_samples

    lmbd_grid_length=int(100*np.max(lmbd_samples))+1
    lmbd_max_rdd=0.01*(lmbd_grid_length-1)
    phi_grid_length=int(100*np.max(phi_samples))+1
    phi_max_rdd=0.01*(phi_grid_length-1)
    N_inv_grid_length=int(100*np.max(N_inv_samples))+1
    N_inv_max_rdd=0.01*(N_inv_grid_length-1)

    H,edges=np.histogramdd(sample_array,
            bins=(np.linspace(0,lmbd_max_rdd,lmbd_grid_length),
                  np.linspace(0,phi_max_rdd,phi_grid_length),
                  np.linspace(0,N_inv_max_rdd,N_inv_grid_length)))
    H_normed=H/np.sum(H)
    mle_loc=np.unravel_index(np.argmax(H_normed),H_normed.shape)
    current_max=H_normed[mle_loc]
    H_normed[mle_loc]==0
    interval_weight=current_max
    print('Now calculating 95% confidence interval.')
    start_time=time.time()
    next_marker=0.1
    while interval_weight<0.95:
        max_loc=np.unravel_index(np.argmax(H_normed),H_normed.shape)
        current_max=H_normed[max_loc]
        H_normed[max_loc]=0
        interval_weight+=current_max
        if interval_weight>next_marker:
            next_marker=np.around(interval_weight,decimals=1)
            if next_marker<interval_weight:
                next_marker+=0.1
            print('Surpassed',next_marker-0.1,'at time',time.time()-start_time,'.')
    H_normed=H/np.sum(H)
    lmbdgrid,phigrid,N_invgrid=np.meshgrid(edges[0],edges[1],edges[2])
    lmbd_ci=[np.min(lmbdgrid[np.where(H_normed>=current_max)]),np.max(lmbdgrid[np.where(H_normed>=current_max)])]
    phi_ci=[np.min(phigrid[np.where(H_normed>=current_max)]),np.max(phigrid[np.where(H_normed>=current_max)])]
    N_inv_ci=[np.min(N_invgrid[np.where(H_normed>=current_max)]),np.max(N_invgrid[np.where(H_normed>=current_max)])]

    return lmbd_mle,lmbd_ci,lmbd_samples,phi_mle,phi_ci,phi_samples,N_inv_mle,N_inv_ci,N_inv_samples
