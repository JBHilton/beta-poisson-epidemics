import numpy as np
from os import mkdir
from os.path import isdir
from pickle import dump
from time import time as get_time
from datasets import cdc_sars_data
from functions import generate_mle_dict, beta_poisson_bootstrap

if isdir('outputs/mles') is False:
    mkdir('outputs/mles')

plague_mles = generate_mle_dict(cdc_sars_data,
                      1.5,
                      0.5,
                      1/2,
                      1.5,
                      0.5)

beta_poi_start=get_time()
beta_poi_lmbd_mle,beta_poi_lmbd_ci,beta_poi_lmbd_samples,beta_poi_phi_mle,beta_poi_phi_ci,beta_poi_phi_samples,beta_poi_N_inv_mle,beta_poi_N_inv_ci,beta_poi_N_inv_samples=beta_poisson_bootstrap(cdc_sars_data,5,1,1/5)
print('Beta-Poisson intervals',beta_poi_lmbd_ci,beta_poi_phi_ci,beta_poi_N_inv_ci,'calculated in',get_time()-beta_poi_start,'seconds.')

bp_mles = plague_mles['beta-Poisson']
print('lambda=', bp_mles[0], beta_poi_lmbd_ci,
      'phi=', bp_mles[1], beta_poi_phi_ci,
      '1/N=', bp_mles[2], beta_poi_N_inv_ci)

# fname = 'outputs/mles/results.pkl'
# with open(fname, 'wb') as f:
#     dump(
#         (plague_mles),
#         f)
