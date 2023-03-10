from argparse import ArgumentParser
from numpy import array, isnan, linspace, mean, meshgrid, percentile, size, var, where, zeros
from os import mkdir
from os.path import isdir, isfile
from pickle import dump, load
from random import choices
from copy import deepcopy
from multiprocessing import Pool
from time import time as get_time
from datasets import (plague_data, mpox_data, nigeria_ebola_data,
    guinea_ebola_data, singapore_sars_data, sk_mers_data, sa_mers_data, noro_data)
from functions import (beta_poisson_loglh, neg_bin_loglh)

data_dict = {
    'plague_data' : plague_data,
    'mpox_data' : mpox_data,
    'nigeria_ebola_data' : nigeria_ebola_data,
    'guinea_ebola_data' : guinea_ebola_data,
    'singapore_sars_data' : singapore_sars_data,
    'sk_mers_data' : sk_mers_data,
    'sa_mers_data' : sa_mers_data,
    'noro_data' :noro_data
    }

if isdir('outputs/sensitivity_analyses') is False:
    mkdir('outputs/sensitivity_analyses')

class LmbdGridCalculator:
    def __init__(self, data_set, mle_dict):
        sample_var = var(data_set)

        self.data_set = data_set
        self.mle_dict = mle_dict
        self.lmbd_mle = mle_dict['beta-Poisson'][0]

    def __call__(self, p):
        try:
            results = self._sensitivity_calculations(p)
        except ValueError as err:
            print(
                'Exception raised for parameters={0}\n\tException: {1}'.format(
                p, err)
                )
            return 0.0
        return results

    def _sensitivity_calculations(self, p):

        phi_p = p[0]
        nu_p = p[1]

        if nu_p <= 1/self.lmbd_mle:
            lmbd_grid_val = beta_poisson_loglh(self.data_set,
                                            self.lmbd_mle,
                                            phi_p,
                                            nu_p)
            return lmbd_grid_val
        else:
            return 1000

class PhiGridCalculator:
    def __init__(self, data_set, mle_dict):
        sample_var = var(data_set)

        self.data_set = data_set
        self.mle_dict = mle_dict
        self.phi_mle = mle_dict['beta-Poisson'][1]

    def __call__(self, p):
        try:
            results = self._sensitivity_calculations(p)
        except ValueError as err:
            print(
                'Exception raised for parameters={0}\n\tException: {1}'.format(
                p, err)
                )
            return 0.0
        return results

    def _sensitivity_calculations(self, p):

        lmbd_p = p[0]
        nu_p = p[1]

        if nu_p <= 1/lmbd_p:
            phi_grid_val = beta_poisson_loglh(self.data_set,
                                          lmbd_p,
                                          self.phi_mle,
                                          nu_p)
            return phi_grid_val
        else:
            return 1000

class NuGridCalculator:
    def __init__(self, data_set, mle_dict):
        sample_var = var(data_set)

        self.data_set = data_set
        self.mle_dict = mle_dict
        self.nu_mle = mle_dict['beta-Poisson'][2]

    def __call__(self, p):
        try:
            results = self._sensitivity_calculations(p)
        except ValueError as err:
            print(
                'Exception raised for parameters={0}\n\tException: {1}'.format(
                p, err)
                )
            return 0.0
        return results

    def _sensitivity_calculations(self, p):

        lmbd_p = p[0]
        phi_p = p[1]

        if self.nu_mle <= 1/lmbd_p:
            nu_grid_val = beta_poisson_loglh(self.data_set,
                                         lmbd_p,
                                         phi_p,
                                         self.nu_mle)
            return nu_grid_val
        else:
            return 1000

def main(no_of_workers,
         data_name):
    main_start = get_time()

    fname = 'outputs/mles/'+data_name+'_results.pkl'
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

    data_set = data_dict[data_name]

    lmbd_results = []
    phi_results = []
    nu_results = []
    lmbd_calculator = LmbdGridCalculator(data_set, mle_dict)
    phi_calculator = PhiGridCalculator(data_set, mle_dict)
    nu_calculator = NuGridCalculator(data_set, mle_dict)

    lmbd_vals = linspace(.05, 5.05, 100)
    phi_vals = linspace(.025, 5.025, 100)
    nu_vals = linspace(0., 1., 100)
    no_vals = len(lmbd_vals) * len(phi_vals) * len(nu_vals)

    lmbd_grid_phi, lmbd_grid_nu = meshgrid(phi_vals, nu_vals)
    phi_grid_lmbd, phi_grid_nu = meshgrid(lmbd_vals, nu_vals)
    nu_grid_lmbd, nu_grid_phi = meshgrid(lmbd_vals, phi_vals)

    no_lmbd_vals = len(phi_vals) * len(nu_vals)
    no_phi_vals = len(lmbd_vals) * len(nu_vals)
    no_nu_vals = len(lmbd_vals) * len(phi_vals)

    lmbd_grid_phi = lmbd_grid_phi.reshape(no_lmbd_vals)
    lmbd_grid_nu = lmbd_grid_nu.reshape(no_lmbd_vals)
    phi_grid_lmbd = phi_grid_lmbd.reshape(no_phi_vals)
    phi_grid_nu = phi_grid_nu.reshape(no_phi_vals)
    nu_grid_lmbd = nu_grid_lmbd.reshape(no_nu_vals)
    nu_grid_phi = nu_grid_phi.reshape(no_nu_vals)

    lmbd_params = [[lmbd_grid_phi[i], lmbd_grid_nu[i]] for i in range(no_lmbd_vals)]
    phi_params = [[phi_grid_lmbd[i], phi_grid_nu[i]] for i in range(no_phi_vals)]
    nu_params = [[nu_grid_lmbd[i], nu_grid_phi[i]] for i in range(no_nu_vals)]

    with Pool(no_of_workers) as pool:
        lmbd_results = pool.map(lmbd_calculator, lmbd_params)
        phi_results = pool.map(phi_calculator, phi_params)
        nu_results = pool.map(nu_calculator, nu_params)


    lmbd_curve = [beta_poisson_loglh(
                    data_set,
                    lmbd_p,
                    mle_dict['beta-Poisson'][1],
                    mle_dict['beta-Poisson'][2]) for lmbd_p in lmbd_vals]
    phi_curve = [beta_poisson_loglh(
                    data_set,
                    mle_dict['beta-Poisson'][0],
                    phi_p,
                    mle_dict['beta-Poisson'][2]) for phi_p in phi_vals]
    nu_curve = [beta_poisson_loglh(
                    data_set,
                    mle_dict['beta-Poisson'][0],
                    mle_dict['beta-Poisson'][1],
                    nu_p) for nu_p in nu_vals]
    lmbd_grid = array([r for r in lmbd_results]).reshape(len(nu_vals), len(phi_vals))
    phi_grid = array([r for r in phi_results]).reshape(len(nu_vals), len(lmbd_vals))
    nu_grid = array([r for r in nu_results]).reshape(len(phi_vals), len(lmbd_vals))

    print('Sensitivity analysis took',get_time()-main_start,'seconds.')

    fname = 'outputs/sensitivity_analyses/'+data_name+'_results.pkl'
    with open(fname, 'wb') as f:
        dump(
            (lmbd_vals,
            phi_vals,
            nu_vals,
            lmbd_curve,
            phi_curve,
            nu_curve,
            lmbd_grid,
            phi_grid,
            nu_grid),
            f)

    return -1

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--no_of_workers', type=int, default=4)
    parser.add_argument('--data_name',
                        type=str,
                        default='plague_data')
    args = parser.parse_args()

    main(args.no_of_workers,
         args.data_name)
