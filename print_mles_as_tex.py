'''This prints the MLE.'''

import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from datasets import (plague_data, monkeypox_data, fasina_ebola_data,
    fay_ebola_data, cdc_sars_data, cowling_mers_data, mers_data, noro_data)
from functions import beta_poisson_pmf, zip_pmf

data_name_list = [
    'plague_data',
    'monkeypox_data',
    'fasina_ebola_data',
    'fay_ebola_data',
    'cdc_sars_data',
    'cowling_mers_data',
    'mers_data',
    'noro_data'
    ]
data_set_list = [
    plague_data,
    monkeypox_data,
    fasina_ebola_data,
    fay_ebola_data,
    cdc_sars_data,
    cowling_mers_data,
    mers_data,
    noro_data
    ]

mle_list = []
ci_list = []
llh_list = []
for i, data_name in enumerate(data_name_list):
    fname = 'outputs/mles/'+data_name+'_results.pkl'
    with open(fname,'rb') as f:
        (mle_dict,
        ci_dict,
        llh_dict) = load(f)
    mle_list.append(mle_dict)
    ci_list.append(ci_dict)
    llh_list.append(llh_dict)

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\toprule
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\
    \midrule
		Plague & ''',mle_list[0]['poisson'],tuple(ci_list[0]['poisson']),''' & ''',mle_list[0]['geometric'],tuple(ci_list[0]['geometric']),''' & ''',mle_list[0]['negative binomial'][0],tuple(ci_list[0]['negative binomial'][0]),''' & ''',mle_list[0]['zip'][0],tuple(ci_list[0]['zip'][0]),''' & ''',mle_list[0]['beta-Poisson'][0],tuple(ci_list[0]['beta-Poisson'][0]),''' \\
		Monkepox & ''',mle_list[1]['poisson'],tuple(ci_list[1]['poisson']),''' & ''',mle_list[1]['geometric'],tuple(ci_list[1]['geometric']),''' & ''',mle_list[1]['negative binomial'][0],tuple(ci_list[1]['negative binomial'][0]),''' & ''',mle_list[1]['zip'][0],tuple(ci_list[1]['zip'][0]),''' & ''',mle_list[1]['beta-Poisson'][0],tuple(ci_list[1]['beta-Poisson'][0]),''' \\
		Ebola, Nigeria 2014 & ''',mle_list[2]['poisson'],tuple(ci_list[2]['poisson']),''' & ''',mle_list[2]['geometric'],tuple(ci_list[2]['geometric']),''' & ''',mle_list[2]['negative binomial'][0],tuple(ci_list[2]['negative binomial'][0]),''' & ''',mle_list[2]['zip'][0],tuple(ci_list[2]['zip'][0]),''' & ''',mle_list[2]['beta-Poisson'][0],tuple(ci_list[2]['beta-Poisson'][0]),''' \\
		Ebola, Guinea 2014 & ''',mle_list[3]['poisson'],tuple(ci_list[3]['poisson']),''' & ''',mle_list[3]['geometric'],tuple(ci_list[3]['geometric']),''' & ''',mle_list[3]['negative binomial'][0],tuple(ci_list[3]['negative binomial'][0]),''' & ''',mle_list[3]['zip'][0],tuple(ci_list[3]['zip'][0]),''' & ''',mle_list[3]['beta-Poisson'][0],tuple(ci_list[3]['beta-Poisson'][0]),''' \\
		SARS, Singapore 2003 & ''',mle_list[4]['poisson'],tuple(ci_list[4]['poisson']),''' & ''',mle_list[4]['geometric'],tuple(ci_list[4]['geometric']),''' & ''',mle_list[4]['negative binomial'][0],tuple(ci_list[4]['negative binomial'][0]),''' & ''',mle_list[4]['zip'][0],tuple(ci_list[4]['zip'][0]),''' & ''',mle_list[4]['beta-Poisson'][0],tuple(ci_list[4]['beta-Poisson'][0]),''' \\
		MERS, South Korea 2015 & ''',mle_list[5]['poisson'],tuple(ci_list[5]['poisson']),''' & ''',mle_list[5]['geometric'],tuple(ci_list[5]['geometric']),''' & ''',mle_list[5]['negative binomial'][0],tuple(ci_list[5]['negative binomial'][0]),''' & ''',mle_list[5]['zip'][0],tuple(ci_list[5]['zip'][0]),''' & ''',mle_list[5]['beta-Poisson'][0],tuple(ci_list[5]['beta-Poisson'][0]),''' \\
		MERS, Saudi Arabia 2015 & ''',mle_list[6]['poisson'],tuple(ci_list[6]['poisson']),''' & ''',mle_list[6]['geometric'],tuple(ci_list[6]['geometric']),''' & ''',mle_list[6]['negative binomial'][0],tuple(ci_list[6]['negative binomial'][0]),''' & ''',mle_list[6]['zip'][0],tuple(ci_list[6]['zip'][0]),''' & ''',mle_list[6]['beta-Poisson'][0],tuple(ci_list[6]['beta-Poisson'][0]),''' \\
		Norovirus, Netherlands 2012 & ''',mle_list[7]['poisson'],tuple(ci_list[7]['poisson']),''' & ''',mle_list[7]['geometric'],tuple(ci_list[7]['geometric']),''' & ''',mle_list[7]['negative binomial'][0],tuple(ci_list[7]['negative binomial'][0]),''' & ''',mle_list[7]['zip'][0],tuple(ci_list[7]['zip'][0]),''' & ''',mle_list[7]['beta-Poisson'][0],tuple(ci_list[7]['beta-Poisson'][0]),''' \\
    \bottomrule
	\end{tabular}
	\caption{Maximum likelihood estimates and confidence intervals (in parentheses) for $\lambda$ by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lccccc}
		\\toprule
		\\textbf{Dataset} & \\textbf{Neg. bin. } $\\theta$ & \\textbf{ZIP} $\sigma$ & \\textbf{Beta-Poisson} $\Phi$ & \\textbf{Beta-Poisson} $N$ \\
    \midrule
		Plague & ''',mle_list[0]['negative binomial'][1],tuple(ci_list[0]['negative binomial'])[1],''' & ''',mle_list[0]['zip'][1],tuple(ci_list[0]['zip'])[1],''' & ''',mle_list[0]['beta-Poisson'][1],tuple(ci_list[0]['beta-Poisson'])[1],''' & ''',mle_list[0]['beta-Poisson'][2],tuple(ci_list[0]['beta-Poisson'])[2],'''\\
		Monkepox & ''',mle_list[1]['negative binomial'][1],tuple(ci_list[1]['negative binomial'])[1],''' & ''',mle_list[1]['zip'][1],tuple(ci_list[1]['zip'])[1],''' & ''',mle_list[1]['beta-Poisson'][1],tuple(ci_list[1]['beta-Poisson'])[1],''' & ''',mle_list[1]['beta-Poisson'][2],tuple(ci_list[1]['beta-Poisson'])[2],'''\\
		Ebola, Nigeria 2014 & ''',mle_list[2]['negative binomial'][1],tuple(ci_list[2]['negative binomial'])[1],''' & ''',mle_list[2]['zip'][1],tuple(ci_list[2]['zip'])[1],''' & ''',mle_list[2]['beta-Poisson'][1],tuple(ci_list[2]['beta-Poisson'])[1],''' & ''',mle_list[2]['beta-Poisson'][2],tuple(ci_list[2]['beta-Poisson'])[2],'''\\
		Ebola, Guinea 2014 & ''',mle_list[3]['negative binomial'][1],tuple(ci_list[3]['negative binomial'])[1],''' & ''',mle_list[3]['zip'][1],tuple(ci_list[3]['zip'])[1],''' & ''',mle_list[3]['beta-Poisson'][1],tuple(ci_list[3]['beta-Poisson'])[1],''' & ''',mle_list[3]['beta-Poisson'][2],tuple(ci_list[3]['beta-Poisson'])[2],'''\\
		SARS, Singapore 2003 & ''',mle_list[4]['negative binomial'][1],tuple(ci_list[4]['negative binomial'])[1],''' & ''',mle_list[4]['zip'][1],tuple(ci_list[4]['zip'])[1],''' & ''',mle_list[4]['beta-Poisson'][1],tuple(ci_list[4]['beta-Poisson'])[1],''' & ''',mle_list[4]['beta-Poisson'][2],tuple(ci_list[4]['beta-Poisson'])[2],'''\\
		MERS, South Korea 2015 & ''',mle_list[5]['negative binomial'][1],tuple(ci_list[5]['negative binomial'])[1],''' & ''',mle_list[5]['zip'][1],tuple(ci_list[5]['zip'])[1],''' & ''',mle_list[5]['beta-Poisson'][1],tuple(ci_list[5]['beta-Poisson'])[1],''' & ''',mle_list[5]['beta-Poisson'][2],tuple(ci_list[5]['beta-Poisson'])[2],'''\\
		MERS, Saudi Arabia 2015 & ''',mle_list[6]['negative binomial'][1],tuple(ci_list[6]['negative binomial'])[1],''' & ''',mle_list[6]['zip'][1],tuple(ci_list[6]['zip'])[1],''' & ''',mle_list[6]['beta-Poisson'][1],tuple(ci_list[6]['beta-Poisson'])[1],''' & ''',mle_list[6]['beta-Poisson'][2],tuple(ci_list[6]['beta-Poisson'])[2],'''\\
		Norovirus, Netherlands 2012 & ''',mle_list[7]['negative binomial'][1],tuple(ci_list[7]['negative binomial'])[1],''' & ''',mle_list[7]['zip'][1],tuple(ci_list[7]['zip'])[1],''' & ''',mle_list[7]['beta-Poisson'][1],tuple(ci_list[7]['beta-Poisson'])[1],''' & ''',mle_list[7]['beta-Poisson'][2],tuple(ci_list[7]['beta-Poisson'])[2],'''\\
	\bottomrule
	\end{tabular}
	\caption{Maximum likelihood estimates and confidence intervals (in parentheses) for model parameters given each dataset.}
\end{table}
''')



print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\toprule
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\
    \midrule
		Plague & ''',llh_list[0]['poisson'][0],''' & ''',llh_list[0]['geometric'][0],''' & ''',llh_list[0]['negative binomial'][0],''' & ''',llh_list[0]['zip'][0],''' & ''',llh_list[0]['beta-Poisson'][0],''' \\
		Monkepox & ''',llh_list[1]['poisson'][0],''' & ''',llh_list[1]['geometric'][0],''' & ''',llh_list[1]['negative binomial'][0],''' & ''',llh_list[1]['zip'][0],''' & ''',llh_list[1]['beta-Poisson'][0],''' \\
		Ebola, Nigeria 2014 & ''',llh_list[2]['poisson'][0],''' & ''',llh_list[2]['geometric'][0],''' & ''',llh_list[2]['negative binomial'][0],''' & ''',llh_list[2]['zip'][0],''' & ''',llh_list[2]['beta-Poisson'][0],''' \\
		Ebola, Guinea 2014 & ''',llh_list[3]['poisson'][0],''' & ''',llh_list[3]['geometric'][0],''' & ''',llh_list[3]['negative binomial'][0],''' & ''',llh_list[3]['zip'][0],''' & ''',llh_list[3]['beta-Poisson'][0],''' \\
		SARS, Singapore 2003 & ''',llh_list[4]['poisson'][0],''' & ''',llh_list[4]['geometric'][0],''' & ''',llh_list[4]['negative binomial'][0],''' & ''',llh_list[4]['zip'][0],''' & ''',llh_list[4]['beta-Poisson'][0],''' \\
		MERS, South Korea 2015 & ''',llh_list[5]['poisson'][0],''' & ''',llh_list[5]['geometric'][0],''' & ''',llh_list[5]['negative binomial'][0],''' & ''',llh_list[5]['zip'][0],''' & ''',llh_list[5]['beta-Poisson'][0],''' \\
		MERS, Saudi Arabia 2015 & ''',llh_list[6]['poisson'][0],''' & ''',llh_list[6]['geometric'][0],''' & ''',llh_list[6]['negative binomial'][0],''' & ''',llh_list[6]['zip'][0],''' & ''',llh_list[6]['beta-Poisson'][0],''' \\
		Norovirus, Netherlands 2012 & ''',llh_list[7]['poisson'][0],''' & ''',llh_list[7]['geometric'][0],''' & ''',llh_list[7]['negative binomial'][0],''' & ''',llh_list[7]['zip'][0],''' & ''',llh_list[7]['beta-Poisson'][0],''' \\
    \bottomrule
	\end{tabular}
	\caption{Log likelihood at MLE by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\toprule
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\
    \midrule
		Plague & ''',llh_list[0]['poisson'][1],''' & ''',llh_list[0]['geometric'][1],''' & ''',llh_list[0]['negative binomial'][1],''' & ''',llh_list[0]['zip'][1],''' & ''',llh_list[0]['beta-Poisson'][1],''' \\
		Monkepox & ''',llh_list[1]['poisson'][1],''' & ''',llh_list[1]['geometric'][1],''' & ''',llh_list[1]['negative binomial'][1],''' & ''',llh_list[1]['zip'][1],''' & ''',llh_list[1]['beta-Poisson'][1],''' \\
		Ebola, Nigeria 2014 & ''',llh_list[2]['poisson'][1],''' & ''',llh_list[2]['geometric'][1],''' & ''',llh_list[2]['negative binomial'][1],''' & ''',llh_list[2]['zip'][1],''' & ''',llh_list[2]['beta-Poisson'][1],''' \\
		Ebola, Guinea 2014 & ''',llh_list[3]['poisson'][1],''' & ''',llh_list[3]['geometric'][1],''' & ''',llh_list[3]['negative binomial'][1],''' & ''',llh_list[3]['zip'][1],''' & ''',llh_list[3]['beta-Poisson'][1],''' \\
		SARS, Singapore 2003 & ''',llh_list[4]['poisson'][1],''' & ''',llh_list[4]['geometric'][1],''' & ''',llh_list[4]['negative binomial'][1],''' & ''',llh_list[4]['zip'][1],''' & ''',llh_list[4]['beta-Poisson'][1],''' \\
		MERS, South Korea 2015 & ''',llh_list[5]['poisson'][1],''' & ''',llh_list[5]['geometric'][1],''' & ''',llh_list[5]['negative binomial'][1],''' & ''',llh_list[5]['zip'][1],''' & ''',llh_list[5]['beta-Poisson'][1],''' \\
		MERS, Saudi Arabia 2015 & ''',llh_list[6]['poisson'][1],''' & ''',llh_list[6]['geometric'][1],''' & ''',llh_list[6]['negative binomial'][1],''' & ''',llh_list[6]['zip'][1],''' & ''',llh_list[6]['beta-Poisson'][1],''' \\
		Norovirus, Netherlands 2012 & ''',llh_list[7]['poisson'][1],''' & ''',llh_list[7]['geometric'][1],''' & ''',llh_list[7]['negative binomial'][1],''' & ''',llh_list[7]['zip'][1],''' & ''',llh_list[7]['beta-Poisson'][1],''' \\
    \bottomrule
	\end{tabular}
	\caption{Akaike information criterion at MLE by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\toprule
		\\textbf{Dataset} & \\textbf{Neg. bin. } & \\textbf{ZIP} \\
    \midrule
		Plague & ''',np.exp(llh_list[0]['negative binomial'][0]-llh_list[0]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[0]['zip'][0]-llh_list[0]['beta-Poisson'][0]),''' \\
		Monkepox & ''',np.exp(llh_list[1]['negative binomial'][0]-llh_list[1]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[1]['zip'][0]-llh_list[1]['beta-Poisson'][0]),''' \\
		Ebola, Nigeria 2014 & ''',np.exp(llh_list[2]['negative binomial'][0]-llh_list[2]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[2]['zip'][0]-llh_list[2]['beta-Poisson'][0]),''' \\
		Ebola, Guinea 2014 & ''',np.exp(llh_list[3]['negative binomial'][0]-llh_list[3]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[3]['zip'][0]-llh_list[3]['beta-Poisson'][0]),''' \\
		SARS, Singapore 2003 & ''',np.exp(llh_list[4]['negative binomial'][0]-llh_list[4]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[4]['zip'][0]-llh_list[4]['beta-Poisson'][0]),''' \\
		MERS, South Korea 2015 & ''',np.exp(llh_list[5]['negative binomial'][0]-llh_list[5]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[5]['zip'][0]-llh_list[5]['beta-Poisson'][0]),''' \\
		MERS, Saudi Arabia 2015 & ''',np.exp(llh_list[6]['negative binomial'][0]-llh_list[6]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[6]['zip'][0]-llh_list[6]['beta-Poisson'][0]),''' \\
		Norovirus, Netherlands 2012 & ''',np.exp(llh_list[7]['negative binomial'][0]-llh_list[7]['beta-Poisson'][0]),''' & ''',np.exp(llh_list[7]['zip'][0]-llh_list[7]['beta-Poisson'][0]),''' \\
	\bottomrule
	\end{tabular}
	\caption{Likelihood ratios of beta-Poisson to negative binomial and ZIP models under each dataset.}
\end{table}
''')
