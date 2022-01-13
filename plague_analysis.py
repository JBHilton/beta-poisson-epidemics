import numpy as np
from os import mkdir
from os.path import isdir
from pickle import dump
from time import time as get_time
from datasets import plague_data
from functions import beta_poisson_loglh, beta_poisson_extinction_prob, empirical_loglh, geom_extinction_prob, geo_loglh, get_phi_and_N_mles, get_theta_mle, get_zip_mles, neg_bin_extinction_prob, neg_bin_loglh, neg_bin_loglh_theta, poisson_loglh, poisson_extinction_prob, zip_loglh, zip_extinction_prob

if isdir('outputs/plague') is False:
    mkdir('outputs/plague')

theta_mle=get_theta_mle(plague_data,1.5)
phi_mle,N_inv_mle=get_phi_and_N_mles(plague_data,0.5,1/2)
lmbd_mle,sigma_mle=get_zip_mles(plague_data,1.5,0.5)

poi_sup=poisson_loglh(plague_data,np.mean(plague_data))
geo_sup=geo_loglh(plague_data,np.mean(plague_data))
nbinom_sup = neg_bin_loglh_theta(plague_data,np.mean(plague_data),theta_mle)
zip_sup = zip_loglh(plague_data,lmbd_mle,sigma_mle)
bp_sup = beta_poisson_loglh(plague_data,np.mean(plague_data),phi_mle,1/N_inv_mle)

e_lh=empirical_loglh(plague_data)

poi_ext=poisson_extinction_prob(np.mean(plague_data))
geom_ext=geom_extinction_prob(np.mean(plague_data))
neg_bin_ext=neg_bin_extinction_prob(np.mean(plague_data),theta_mle)
zip_ext=zip_extinction_prob(lmbd_mle,sigma_mle)
beta_poi_ext=beta_poisson_extinction_prob(np.mean(plague_data),phi_mle,1/N_inv_mle)

true_lmbd,true_phi,true_N=np.mean(plague_data),phi_mle,1/N_inv_mle
data=plague_data
lmbd_lh_array=beta_poisson_loglh(data,np.linspace(0.01,5,500),true_phi,true_N)
phi_lh_array=np.zeros(1000)
for i in range(1000):
    phi_lh_array[i]=beta_poisson_loglh(data,true_lmbd,(i+1)*1e-2,true_N)
n_lh_array=np.zeros(750)
n_lh_array[0]=neg_bin_loglh(data,true_lmbd,true_phi)
for i in range(1,750):
    n_lh_array[i]=beta_poisson_loglh(data,true_lmbd,true_phi,1/(i/1000))

data = plague_data
true_lmbd,true_phi,true_N=np.mean(plague_data),phi_mle,1/N_inv_mle

lambda_range = np.linspace(0.1,5,10)
phi_range = np.linspace(0.1,10,10)
nu_range = np.linspace(0,0.75,10)

lambda_fixed_size = len(phi_range)*len(nu_range)
phi_vect,nu_vect = np.meshgrid(phi_range,nu_range)
phi_vect = phi_vect.reshape(lambda_fixed_size,1)
nu_vect = nu_vect.reshape(lambda_fixed_size,1)
loglh_grid_lambda_fixed = np.zeros((lambda_fixed_size,1))
start_time = get_time()
for i in range(lambda_fixed_size):
    if nu_vect[i]==0:
        loglh_grid_lambda_fixed[i] = neg_bin_loglh(data,true_lmbd,phi_vect[i])
    else:
        loglh_grid_lambda_fixed[i] = beta_poisson_loglh(data,true_lmbd,phi_vect[i],1/nu_vect[i])
    if (i%10000==0)&(i>0):
        print(i,'iterations completed in',get_time()-start_time,'seconds,',lambda_fixed_size-i,'iterations left, estimated',(lambda_fixed_size-i)*(get_time()-start_time)/i,'seconds left.')

print('Calculation with lambda fixed completed in',get_time()-start_time,'seconds.')

phi_fixed_size = len(lambda_range)*len(nu_range)
lambda_vect,nu_vect = np.meshgrid(lambda_range,nu_range)
lambda_vect = lambda_vect.reshape(phi_fixed_size,1)
nu_vect = nu_vect.reshape(phi_fixed_size,1)
loglh_grid_phi_fixed = np.zeros((phi_fixed_size,1))
start_time = get_time()
for i in range(phi_fixed_size):
    if nu_vect[i]==0:
        loglh_grid_phi_fixed[i] = neg_bin_loglh(data,lambda_vect[i],true_phi)
    else:
        loglh_grid_phi_fixed[i] = beta_poisson_loglh(data,lambda_vect[i],true_phi,1/nu_vect[i])
    if (i%10000==0)&(i>0):
        print(i,'iterations completed in',get_time()-start_time,'seconds,',phi_fixed_size-i,'iterations left, estimated',(phi_fixed_size-i)*(get_time()-start_time)/i,'seconds left.')

print('Calculation with phi fixed completed in',get_time()-start_time,'seconds.')

nu_fixed_size = len(lambda_range)*len(phi_range)
lambda_vect,phi_vect = np.meshgrid(lambda_range,phi_range)
lambda_vect = lambda_vect.reshape(nu_fixed_size,1)
phi_vect = phi_vect.reshape(nu_fixed_size,1)
loglh_grid_nu_fixed = np.zeros((nu_fixed_size,1))
start_time = get_time()
for i in range(nu_fixed_size):
    loglh_grid_nu_fixed[i] = beta_poisson_loglh(data,lambda_vect[i],phi_vect[i],true_N)
    if (i%10000==0)&(i>0):
        print(i,'iterations completed in',get_time()-start_time,'seconds,',nu_fixed_size-i,'iterations left, estimated',(nu_fixed_size-i)*(get_time()-start_time)/i,'seconds left.')

print('Calculation with nu fixed completed in',get_time()-start_time,'seconds.')

fname = 'outputs/plague/results.pkl'
with open(fname, 'wb') as f:
    dump(
        (theta_mle,
        lmbd_mle,
        phi_mle,
        N_inv_mle,
        sigma_mle),
        f)
