from argparse import ArgumentParser
from numpy import arange, array, mean, meshgrid, percentile, size, var, zeros
from os import mkdir
from os.path import isdir, isfile
from pickle import dump, load
from random import choices
from copy import deepcopy
from multiprocessing import Pool
from time import time as get_time
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)
from functions import (beta_poisson_loglh, neg_bin_loglh)

data_dict = {
    'plague_data' : plague_data,
    'monkeypox_data' : monkeypox_data,
    'fasina_ebola_data' : fasina_ebola_data,
    'fay_ebola_data' : fay_ebola_data,
    'cdc_sars_data' : cdc_sars_data,
    'cowling_mers_data' : cowling_mers_data,
    'mers_data' : mers_data,
    'noro_data' :noro_data
    }

if isdir('outputs/sensitivity_analyses') is False:
    mkdir('outputs/sensitivity_analyses')

class LLHCalculator:
    def __init__(self, data_set, mle_dict):
        sample_var = var(data_set)

        self.data_set = data_set
        self.mle_dict = mle_dict
        self.lmbd_mle = mle_dict['beta-Poisson'][0]
        self.phi_mle = mle_dict['beta-Poisson'][1]
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
        nu_p = p[2]

        lmbd_grid_val = beta_poisson_loglh(self.data_set,
                                           self.lmbd_mle,
                                           phi_p,
                                           nu_p)
        phi_grid_val = beta_poisson_loglh(self.data_set,
                                          lmbd_p,
                                          self.phi_mle,
                                          nu_p)
        nu_grid_val = beta_poisson_loglh(self.data_set,
                                         lmbd_p,
                                         phi_p,
                                         self.nu_mle)

        return [lmbd_grid_val,
                phi_grid_val,
                nu_grid_val,
                p]

def main(no_of_workers,
         data_name):
    main_start = get_time()

    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname,'rb') as f:
        (mle_dict,
        var_dict,
        ci_dict,
        var_ci_dict,
        llh_dict) = load(f)

    data_set = data_dict[data_name]

    results = []
    calculator = LLHCalculator(data_set, mle_dict)

    lmbd_vals = arange(.1, 5., .1)
    phi_vals = arange(.1, 10., .1)
    nu_vals = arange(0., 1., .1)
    no_vals = len(lmbd_vals) * len(phi_vals) * len(nu_vals)

    L, P, Nu = meshgrid(lmbd_vals, phi_vals, nu_vals)
    L = L.reshape(no_vals)
    P = P.reshape(no_vals)
    Nu = Nu.reshape(no_vals)

    params = [[L[i], P[i], Nu[i]] for i in range(no_vals)]

    with Pool(no_of_workers) as pool:
        results = pool.map(calculator, params)


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
    lmbd_grid = [r[0] for r in results]
    phi_grid = [r[1] for r in results]
    nu_grid = [r[2] for r in results]

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
