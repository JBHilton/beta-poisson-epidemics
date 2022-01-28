from argparse import ArgumentParser
from numpy import arange, array, percentile, size
from os import mkdir
from os.path import isdir, isfile
from pickle import dump
from random import choices
from copy import deepcopy
from multiprocessing import Pool
from time import time as get_time
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_sars_data, mers_data, noro_data)
from functions import ci_from_bootstrap_samples, generate_mle_dict

data_dict = {
    'plague_data' : plague_data,
    'monkeypox_data' : monkeypox_data,
    'fasina_ebola_data' : fasina_ebola_data,
    'fay_ebola_data' : fay_ebola_data,
    'cdc_sars_data' : cdc_sars_data,
    'cowling_sars_data' : cowling_sars_data,
    'mers_data' : mers_data,
    'noro_data' :noro_data
    }

confidence_level = 95

if isdir('outputs/mles') is False:
    mkdir('outputs/mles')

class MLECalculator:
    def __init__(self, data_set):
        self.data_set = data_set
        self.mle_dict = generate_mle_dict(data_set,
                              1.5,
                              0.5,
                              1/2,
                              1.5,
                              0.5)
        self.sample_size = sample_size=size(data_set)

    def __call__(self, p):
        try:
            result = self._resample_and_fit(p)
        except ValueError as err:
            print(
                'Exception raised for parameters={0}\n\tException: {1}'.format(
                p, err)
                )
            return 0.0
        return result

    def _resample_and_fit(self, p):
        data_now=choices(self.data_set, k=self.sample_size)
        sample_dict = generate_mle_dict(data_now,
                              1.5,
                              0.5,
                              1/2,
                              1.5,
                              0.5)
        return sample_dict

def main(no_of_workers,
         no_samples,
         data_name):
    main_start = get_time()

    data_set = data_dict[data_name]

    results = []
    calculator = MLECalculator(data_set)
    params = arange(no_samples)

    with Pool(no_of_workers) as pool:
        dict_samples = pool.map(calculator, params)

    mle_dict = calculator.mle_dict

    poisson_lmbd_samples = array([d['poisson'] for d in dict_samples])
    poisson_ci = ci_from_bootstrap_samples(poisson_lmbd_samples,
                                           mle_dict['poisson'],
                                           confidence_level)

    geo_lmbd_samples = array([d['geometric'] for d in dict_samples])
    geo_ci = ci_from_bootstrap_samples(geo_lmbd_samples,
                                       mle_dict['geometric'],
                                       confidence_level)

    neg_bin_lmbd_samples = array(
                            [d['negative binomial'][0] for d in dict_samples])
    neg_bin_lmbd_ci = ci_from_bootstrap_samples(neg_bin_lmbd_samples,
                                            mle_dict['negative binomial'][0],
                                            confidence_level)
    neg_bin_theta_samples = array(
                            [d['negative binomial'][1] for d in dict_samples])
    neg_bin_theta_ci = ci_from_bootstrap_samples(neg_bin_theta_samples,
                                            mle_dict['negative binomial'][1],
                                            confidence_level)

    zip_lmbd_samples = array(
                            [d['zip'][0] for d in dict_samples])
    zip_lmbd_ci = ci_from_bootstrap_samples(zip_lmbd_samples,
                                            mle_dict['zip'][0],
                                            confidence_level)
    zip_sigma_samples = array(
                            [d['zip'][1] for d in dict_samples])
    zip_sigma_ci = ci_from_bootstrap_samples(zip_sigma_samples,
                                            mle_dict['zip'][1],
                                            confidence_level)

    beta_poi_lmbd_samples = array(
                            [d['beta-Poisson'][0] for d in dict_samples])
    beta_poi_lmbd_ci = ci_from_bootstrap_samples(beta_poi_lmbd_samples,
                                                mle_dict['beta-Poisson'][0],
                                                confidence_level)
    beta_poi_phi_samples = array(
                            [d['beta-Poisson'][1] for d in dict_samples])
    beta_poi_phi_ci = ci_from_bootstrap_samples(beta_poi_phi_samples,
                                                mle_dict['beta-Poisson'][1],
                                                confidence_level)
    beta_poi_N_inv_samples = array(
                            [d['beta-Poisson'][2] for d in dict_samples])
    beta_poi_N_inv_ci = ci_from_bootstrap_samples(beta_poi_N_inv_samples,
                                                  mle_dict['beta-Poisson'][2],
                                                  confidence_level)

    ci_dict = {
        'poisson' : poisson_ci,
        'geometric' : geo_ci,
        'negative binomial' : [neg_bin_lmbd_ci, neg_bin_theta_ci],
        'zip' : [zip_lmbd_ci, zip_sigma_ci],
        'beta-Poisson' : [beta_poi_lmbd_ci,
                          beta_poi_phi_ci,
                          beta_poi_N_inv_ci]
    }

    print(calculator.mle_dict)
    print(ci_dict)


    print('Bootstrap took',get_time()-main_start,'seconds.')

    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname, 'wb') as f:
        dump(
            (calculator.mle_dict,
            ci_dict),
            f)

    return -1

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--no_of_workers', type=int, default=8)
    parser.add_argument('--no_samples',
                        type=int,
                        default=1000)
    parser.add_argument('--data_name',
                        type=str,
                        default='plague_data')
    args = parser.parse_args()

    main(args.no_of_workers,
         args.no_samples,
         args.data_name)
