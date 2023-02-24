from argparse import ArgumentParser
from numpy import arange, array, mean, percentile, size, var
from os import mkdir
from os.path import isdir, isfile
from pickle import dump
from random import choices
from copy import deepcopy
from multiprocessing import Pool
from time import time as get_time
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)
from functions import (ci_from_bootstrap_samples, generate_llh_dict,
    generate_mle_dict, generate_p0_dict, generate_superspread_dict, generate_var_dict)

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

confidence_level = 95

if isdir('outputs/mles') is False:
    mkdir('outputs/mles')

class MLECalculator:
    def __init__(self, data_set):

        sample_mean = mean(data_set)
        sample_var = var(data_set)

        # Define some ansatzes of parameter MLEs
        theta_0 = (sample_var / sample_mean) - 1
        phi_0 = 1 / theta_0
        N_0 = 10 * sample_mean
        zip_lmbd_0 = sample_mean
        sigma_0 = 0.5

        self.data_set = data_set
        self.mle_dict = generate_mle_dict(data_set,
                            theta_0,
                            phi_0,
                            1 / N_0,
                            zip_lmbd_0,
                            sigma_0)
        self.var_dict = generate_var_dict(data_set,
                            self.mle_dict)
        self.superspread_dict = generate_superspread_dict(data_set,
                            self.mle_dict)
        self.p0_dict = generate_p0_dict(data_set,
                            self.mle_dict)
        self.sample_size = size(data_set)

    def __call__(self, p):
        try:
            results = self._resample_and_fit(p)
        except ValueError as err:
            print(
                'Exception raised for parameters={0}\n\tException: {1}'.format(
                p, err)
                )
            return 0.0
        return results

    def _resample_and_fit(self, p):
        sample_flag = 0
        while sample_flag==0:
            data_now=choices(self.data_set, k=self.sample_size)
            sample_mean = mean(data_now)
            if sample_mean > 0:
                sample_flag = 1
        sample_var = var(data_now)
        if sample_var > sample_mean:
            theta_0 = (sample_var / sample_mean) - 1
        else:
            theta_0 = 1e-1
        phi_0 = 1 / theta_0
        N_0 = 2 * sample_mean
        zip_lmbd_0 = sample_mean
        sigma_0 = 0.5
        sample_dict = generate_mle_dict(data_now,
                              theta_0,
                              phi_0,
                              1 / N_0,
                              zip_lmbd_0,
                              sigma_0)
        sample_vars = generate_var_dict(data_now, sample_dict)
        sample_props = generate_superspread_dict(data_now,
                            sample_dict)
        sample_p0 = generate_p0_dict(data_now, sample_dict)
        return [sample_dict, sample_vars, sample_props, sample_p0]

def main(no_of_workers,
         no_samples,
         data_name):
    main_start = get_time()

    data_set = data_dict[data_name]

    results = []
    calculator = MLECalculator(data_set)
    params = arange(no_samples)

    with Pool(no_of_workers) as pool:
        results = pool.map(calculator, params)

    dict_samples = [r[0] for r in results]
    var_samples = [r[1] for r in results]
    superspread_samples = [r[2] for r in results]
    p0_samples = [r[3] for r in results]
    mle_dict = calculator.mle_dict
    var_dict = calculator.var_dict
    od_dict = {
        'sample' : (var(data_set) - mean(data_set))/mean(data_set),
        'geometric' : (var_dict['geometric'] - mle_dict['geometric'])/ mle_dict['geometric'],
        'negative binomial' : mle_dict['negative binomial'][1],
        'zip' : (var_dict['zip'] - mle_dict['zip'][0]*(1-mle_dict['zip'][1]))/ (mle_dict['zip'][0] * (1-mle_dict['zip'][1])),
        'beta-Poisson' : (var_dict['beta-Poisson'] - mle_dict['beta-Poisson'][0])/ mle_dict['beta-Poisson'][0]
    }
    superspread_dict = calculator.superspread_dict
    p0_dict = calculator.p0_dict

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
    beta_poi_nu_samples = array(
                            [d['beta-Poisson'][2] for d in dict_samples])
    beta_poi_nu_ci = ci_from_bootstrap_samples(beta_poi_nu_samples,
                                                  mle_dict['beta-Poisson'][2],
                                                  confidence_level)

    ci_dict = {
        'poisson' : poisson_ci,
        'geometric' : geo_ci,
        'negative binomial' : [neg_bin_lmbd_ci, neg_bin_theta_ci],
        'zip' : [zip_lmbd_ci, zip_sigma_ci],
        'beta-Poisson' : [beta_poi_lmbd_ci,
                          beta_poi_phi_ci,
                          beta_poi_nu_ci]
    }

    poisson_var_samples = array([d['poisson'] for d in var_samples])
    poisson_var_ci = ci_from_bootstrap_samples(poisson_var_samples,
                                           var_dict['poisson'],
                                           confidence_level)
    geo_var_samples = array([d['geometric'] for d in var_samples])
    geo_var_ci = ci_from_bootstrap_samples(geo_var_samples,
                                           var_dict['geometric'],
                                           confidence_level)
    neg_bin_var_samples = array([d['negative binomial'] for d in var_samples])
    neg_bin_var_ci = ci_from_bootstrap_samples(neg_bin_var_samples,
                                           var_dict['negative binomial'],
                                           confidence_level)
    zip_var_samples = array([d['zip'] for d in var_samples])
    zip_var_ci = ci_from_bootstrap_samples(zip_var_samples,
                                           var_dict['zip'],
                                           confidence_level)
    beta_poi_var_samples = array([d['beta-Poisson'] for d in var_samples])
    beta_poi_var_ci = ci_from_bootstrap_samples(beta_poi_var_samples,
                                           var_dict['beta-Poisson'],
                                           confidence_level)

    var_ci_dict = {
        'poisson' : poisson_var_ci,
        'geometric' : geo_var_ci,
        'negative binomial' : neg_bin_var_ci,
        'zip' : zip_var_ci,
        'beta-Poisson' : beta_poi_var_ci
    }

    geo_od_samples = (array([d['geometric'] for d in var_samples]) - array([d['geometric'] for d in dict_samples]))/ array([d['geometric'] for d in dict_samples])
    geo_od_ci = ci_from_bootstrap_samples(geo_od_samples,
                                           0,
                                           confidence_level)
    neg_bin_od_samples = array([d['negative binomial'][1] for d in dict_samples])
    neg_bin_od_ci = ci_from_bootstrap_samples(neg_bin_od_samples,
                                           0,
                                           confidence_level)
    zip_od_samples = (array([d['zip'] for d in var_samples]) - array([d['zip'][0]*(1-d['zip'][1]) for d in dict_samples]))/ array([d['zip'][0]*(1-d['zip'][1]) for d in dict_samples])
    zip_od_ci = ci_from_bootstrap_samples(zip_od_samples,
                                           0,
                                           confidence_level)
    beta_poi_od_samples = (array([d['beta-Poisson'] for d in var_samples]) - array([d['beta-Poisson'][0] for d in dict_samples]))/ array([d['beta-Poisson'][0] for d in dict_samples])
    beta_poi_od_ci = ci_from_bootstrap_samples(beta_poi_od_samples,
                                           0,
                                           confidence_level)

    od_ci_dict = {
        'geometric' : geo_od_ci,
        'negative binomial' : neg_bin_od_ci,
        'zip' : zip_od_ci,
        'beta-Poisson' : beta_poi_od_ci
    }

    poisson_superspread_samples = array([d['poisson'] for d in superspread_samples])
    poisson_superspread_ci = ci_from_bootstrap_samples(poisson_superspread_samples,
                                           superspread_dict['poisson'],
                                           confidence_level)
    geo_superspread_samples = array([d['geometric'] for d in superspread_samples])
    geo_superspread_ci = ci_from_bootstrap_samples(geo_superspread_samples,
                                           superspread_dict['geometric'],
                                           confidence_level)
    neg_bin_superspread_samples = array([d['negative binomial'] for d in superspread_samples])
    neg_bin_superspread_ci = ci_from_bootstrap_samples(neg_bin_superspread_samples,
                                           superspread_dict['negative binomial'],
                                           confidence_level)
    zip_superspread_samples = array([d['zip'] for d in superspread_samples])
    zip_superspread_ci = ci_from_bootstrap_samples(zip_superspread_samples,
                                           superspread_dict['zip'],
                                           confidence_level)
    beta_poi_superspread_samples = array([d['beta-Poisson'] for d in superspread_samples])
    beta_poi_superspread_ci = ci_from_bootstrap_samples(beta_poi_superspread_samples,
                                           superspread_dict['beta-Poisson'],
                                           confidence_level)

    superspread_ci_dict = {
        'poisson' : poisson_superspread_ci,
        'geometric' : geo_superspread_ci,
        'negative binomial' : neg_bin_superspread_ci,
        'zip' : zip_superspread_ci,
        'beta-Poisson' : beta_poi_superspread_ci
    }

    poisson_p0_samples = array([d['poisson'] for d in p0_samples])
    poisson_p0_ci = ci_from_bootstrap_samples(poisson_p0_samples,
                                           p0_dict['poisson'],
                                           confidence_level)
    geo_p0_samples = array([d['geometric'] for d in p0_samples])
    geo_p0_ci = ci_from_bootstrap_samples(geo_p0_samples,
                                           p0_dict['geometric'],
                                           confidence_level)
    neg_bin_p0_samples = array([d['negative binomial'] for d in p0_samples])
    neg_bin_p0_ci = ci_from_bootstrap_samples(neg_bin_p0_samples,
                                           p0_dict['negative binomial'],
                                           confidence_level)
    zip_p0_samples = array([d['zip'] for d in p0_samples])
    zip_p0_ci = ci_from_bootstrap_samples(zip_p0_samples,
                                           p0_dict['zip'],
                                           confidence_level)
    beta_poi_p0_samples = array([d['beta-Poisson'] for d in p0_samples])
    beta_poi_p0_ci = ci_from_bootstrap_samples(beta_poi_p0_samples,
                                           p0_dict['beta-Poisson'],
                                           confidence_level)

    p0_ci_dict = {
        'poisson' : poisson_p0_ci,
        'geometric' : geo_p0_ci,
        'negative binomial' : neg_bin_p0_ci,
        'zip' : zip_p0_ci,
        'beta-Poisson' : beta_poi_p0_ci
    }

    llh_dict = generate_llh_dict(data_set, mle_dict)


    print('Bootstrap took',get_time()-main_start,'seconds.')

    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname, 'wb') as f:
        dump(
            (mle_dict,
            var_dict,
            od_dict,
            superspread_dict,
            p0_dict,
            ci_dict,
            var_ci_dict,
            od_ci_dict,
            superspread_ci_dict,
            p0_ci_dict,
            llh_dict),
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
