'''This prints the outputs of our analysis as Latex tables.'''

import matplotlib.pyplot as plt
import numpy as np
from pickle import load
from scipy import stats
from scipy.stats.distributions import chi2
from datasets import (plague_data, mpox_data, nigeria_ebola_data,
    guinea_ebola_data, singapore_sars_data, sk_mers_data, sa_mers_data, noro_data)
from functions import beta_poisson_pmf, zip_pmf

np.set_printoptions(precision=2)

data_name_list = [
    'plague_data',
    'mpox_data',
    'nigeria_ebola_data',
    'guinea_ebola_data',
    'singapore_sars_data',
    'sk_mers_data',
    'sa_mers_data',
    'noro_data'
    ]
data_set_list = [
    plague_data,
    mpox_data,
    nigeria_ebola_data,
    guinea_ebola_data,
    singapore_sars_data,
    sk_mers_data,
    sa_mers_data,
    noro_data
    ]

mle_list = []
ci_list = []
mean_ci_list = []
var_list = []
var_ci_list = []
superspread_list = []
superspread_ci_list = []
p0_list = []
p0_ci_list = []
od_list = []
od_ci_list = []
llh_list = []
for i, data_name in enumerate(data_name_list):
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
    mle_list.append(mle_dict)
    ci_list.append(ci_dict)
    mean_ci_list.append(mean_ci_dict)
    var_list.append(var_dict)
    var_ci_list.append(var_ci_dict)
    superspread_list.append(superspread_dict)
    superspread_ci_list.append(superspread_ci_dict)
    p0_list.append(p0_dict)
    p0_ci_list.append(p0_ci_dict)
    od_list.append(od_dict)
    od_ci_list.append(od_ci_dict)
    llh_list.append(llh_dict)
    
unique_offspring_vals = np.unique(np.hstack(data_set_list))

table_str = '''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{l|c|c|c|c|c|c|c|c}
		\\hline
		\\ & \\textbf{Plague} & \\textbf{Mpox} & \\textbf{Ebola, Nigeria 2014} & \\textbf{Ebola, Guinea 2014} & \\textbf{SARS, Singapore 2003} & \\textbf{MERS, South Korea 2015} & \\textbf{MERS, Saudi Arabia 2015} & \\textbf{Norovirus, Netherlands 2012} \\\\
'''
for v in unique_offspring_vals:
    table_str += (''' \\hline ''' + str(v))
    for i in range(len(data_set_list)):
        ds = data_set_list[i]
        if v==superspread_list[i]['boundary']:
        	table_str += (''' & \\textbf{''' + str(ds.count(v)) + '''}''')
		else:
            table_str += (''' & ''' + str(ds.count(v)))
    table_str += '''\\\\\n'''
table_str += (''' \\hline Total''')
for ds in data_set_list:
        table_str += (''' & ''' + str(sum(ds)))
table_str += '''\\\\\n'''
table_str += '''
    \\hline
	\end{tabular}
	\caption{Frequency of secondary case numbers by dataset.}
\end{table}
'''

print(table_str)

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['poisson'], 2),tuple(np.round(ci_list[0]['poisson'], 2)),''' & ''',np.round(mle_list[0]['geometric'], 2),tuple(np.round(ci_list[0]['geometric'], 2)),''' & ''',np.round(mle_list[0]['negative binomial'][0], 2),tuple(np.round(ci_list[0]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[0]['zip'][0], 2),tuple(np.round(ci_list[0]['zip'][0], 2)),''' & ''',np.round(mle_list[0]['beta-Poisson'][0], 2),tuple(np.round(ci_list[0]['beta-Poisson'][0], 2)),''' \\\\
		Mpox & ''',np.round(mle_list[1]['poisson'], 2),tuple(np.round(ci_list[1]['poisson'], 2)),''' & ''',np.round(mle_list[1]['geometric'], 2),tuple(np.round(ci_list[1]['geometric'], 2)),''' & ''',np.round(mle_list[1]['negative binomial'][0], 2),tuple(np.round(ci_list[1]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[1]['zip'][0], 2),tuple(np.round(ci_list[1]['zip'][0], 2)),''' & ''',np.round(mle_list[1]['beta-Poisson'][0], 2),tuple(np.round(ci_list[1]['beta-Poisson'][0], 2)),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['poisson'], 2),tuple(np.round(ci_list[2]['poisson'], 2)),''' & ''',np.round(mle_list[2]['geometric'], 2),tuple(np.round(ci_list[2]['geometric'], 2)),''' & ''',np.round(mle_list[2]['negative binomial'][0], 2),tuple(np.round(ci_list[2]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[2]['zip'][0], 2),tuple(np.round(ci_list[2]['zip'][0], 2)),''' & ''',np.round(mle_list[2]['beta-Poisson'][0], 2),tuple(np.round(ci_list[2]['beta-Poisson'][0], 2)),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['poisson'], 2),tuple(np.round(ci_list[3]['poisson'], 2)),''' & ''',np.round(mle_list[3]['geometric'], 2),tuple(np.round(ci_list[3]['geometric'], 2)),''' & ''',np.round(mle_list[3]['negative binomial'][0], 2),tuple(np.round(ci_list[3]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[3]['zip'][0], 2),tuple(np.round(ci_list[3]['zip'][0], 2)),''' & ''',np.round(mle_list[3]['beta-Poisson'][0], 2),tuple(np.round(ci_list[3]['beta-Poisson'][0], 2)),''' \\\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['poisson'], 2),tuple(np.round(ci_list[4]['poisson'], 2)),''' & ''',np.round(mle_list[4]['geometric'], 2),tuple(np.round(ci_list[4]['geometric'], 2)),''' & ''',np.round(mle_list[4]['negative binomial'][0], 2),tuple(np.round(ci_list[4]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[4]['zip'][0], 2),tuple(np.round(ci_list[4]['zip'][0], 2)),''' & ''',np.round(mle_list[4]['beta-Poisson'][0], 2),tuple(np.round(ci_list[4]['beta-Poisson'][0], 2)),''' \\\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['poisson'], 2),tuple(np.round(ci_list[5]['poisson'], 2)),''' & ''',np.round(mle_list[5]['geometric'], 2),tuple(np.round(ci_list[5]['geometric'], 2)),''' & ''',np.round(mle_list[5]['negative binomial'][0], 2),tuple(np.round(ci_list[5]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[5]['zip'][0], 2),tuple(np.round(ci_list[5]['zip'][0], 2)),''' & ''',np.round(mle_list[5]['beta-Poisson'][0], 2),tuple(np.round(ci_list[5]['beta-Poisson'][0], 2)),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['poisson'], 2),tuple(np.round(ci_list[6]['poisson'], 2)),''' & ''',np.round(mle_list[6]['geometric'], 2),tuple(np.round(ci_list[6]['geometric'], 2)),''' & ''',np.round(mle_list[6]['negative binomial'][0], 2),tuple(np.round(ci_list[6]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[6]['zip'][0], 2),tuple(np.round(ci_list[6]['zip'][0], 2)),''' & ''',np.round(mle_list[6]['beta-Poisson'][0], 2),tuple(np.round(ci_list[6]['beta-Poisson'][0], 2)),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['poisson'], 2),tuple(np.round(ci_list[7]['poisson'], 2)),''' & ''',np.round(mle_list[7]['geometric'], 2),tuple(np.round(ci_list[7]['geometric'], 2)),''' & ''',np.round(mle_list[7]['negative binomial'][0], 2),tuple(np.round(ci_list[7]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[7]['zip'][0], 2),tuple(np.round(ci_list[7]['zip'][0], 2)),''' & ''',np.round(mle_list[7]['beta-Poisson'][0], 2),tuple(np.round(ci_list[7]['beta-Poisson'][0], 2)),''' \\\\
    \\hline
	\end{tabular}
	\caption{Maximum likelihood estimates and confidence intervals (in parentheses) for $\lambda$ by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['poisson'], 2),tuple(np.round(ci_list[0]['poisson'], 2)),''' & ''',np.round(mle_list[0]['geometric'], 2),tuple(np.round(ci_list[0]['geometric'], 2)),''' & ''',np.round(mle_list[0]['negative binomial'][0], 2),tuple(np.round(ci_list[0]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[0]['zip'][0]*(1-mle_list[0]['zip'][1]), 2),tuple(np.round(mean_ci_list[0]['zip'], 2)),''' & ''',np.round(mle_list[0]['beta-Poisson'][0], 2),tuple(np.round(ci_list[0]['beta-Poisson'][0], 2)),''' \\\\
		Mpox & ''',np.round(mle_list[1]['poisson'], 2),tuple(np.round(ci_list[1]['poisson'], 2)),''' & ''',np.round(mle_list[1]['geometric'], 2),tuple(np.round(ci_list[1]['geometric'], 2)),''' & ''',np.round(mle_list[1]['negative binomial'][0], 2),tuple(np.round(ci_list[1]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[1]['zip'][0]*(1-mle_list[1]['zip'][1]), 2),tuple(np.round(mean_ci_list[1]['zip'], 2)),''' & ''',np.round(mle_list[1]['beta-Poisson'][0], 2),tuple(np.round(ci_list[1]['beta-Poisson'][0], 2)),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['poisson'], 2),tuple(np.round(ci_list[2]['poisson'], 2)),''' & ''',np.round(mle_list[2]['geometric'], 2),tuple(np.round(ci_list[2]['geometric'], 2)),''' & ''',np.round(mle_list[2]['negative binomial'][0], 2),tuple(np.round(ci_list[2]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[2]['zip'][0]*(1-mle_list[2]['zip'][1]), 2),tuple(np.round(mean_ci_list[2]['zip'], 2)),''' & ''',np.round(mle_list[2]['beta-Poisson'][0], 2),tuple(np.round(ci_list[2]['beta-Poisson'][0], 2)),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['poisson'], 2),tuple(np.round(ci_list[3]['poisson'], 2)),''' & ''',np.round(mle_list[3]['geometric'], 2),tuple(np.round(ci_list[3]['geometric'], 2)),''' & ''',np.round(mle_list[3]['negative binomial'][0], 2),tuple(np.round(ci_list[3]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[3]['zip'][0]*(1-mle_list[3]['zip'][1]), 2),tuple(np.round(mean_ci_list[3]['zip'], 2)),''' & ''',np.round(mle_list[3]['beta-Poisson'][0], 2),tuple(np.round(ci_list[3]['beta-Poisson'][0], 2)),''' \\\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['poisson'], 2),tuple(np.round(ci_list[4]['poisson'], 2)),''' & ''',np.round(mle_list[4]['geometric'], 2),tuple(np.round(ci_list[4]['geometric'], 2)),''' & ''',np.round(mle_list[4]['negative binomial'][0], 2),tuple(np.round(ci_list[4]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[4]['zip'][0]*(1-mle_list[4]['zip'][1]), 2),tuple(np.round(mean_ci_list[4]['zip'], 2)),''' & ''',np.round(mle_list[4]['beta-Poisson'][0], 2),tuple(np.round(ci_list[4]['beta-Poisson'][0], 2)),''' \\\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['poisson'], 2),tuple(np.round(ci_list[5]['poisson'], 2)),''' & ''',np.round(mle_list[5]['geometric'], 2),tuple(np.round(ci_list[5]['geometric'], 2)),''' & ''',np.round(mle_list[5]['negative binomial'][0], 2),tuple(np.round(ci_list[5]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[5]['zip'][0]*(1-mle_list[5]['zip'][1]), 2),tuple(np.round(mean_ci_list[5]['zip'], 2)),''' & ''',np.round(mle_list[5]['beta-Poisson'][0], 2),tuple(np.round(ci_list[5]['beta-Poisson'][0], 2)),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['poisson'], 2),tuple(np.round(ci_list[6]['poisson'], 2)),''' & ''',np.round(mle_list[6]['geometric'], 2),tuple(np.round(ci_list[6]['geometric'], 2)),''' & ''',np.round(mle_list[6]['negative binomial'][0], 2),tuple(np.round(ci_list[6]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[6]['zip'][0]*(1-mle_list[6]['zip'][1]), 2),tuple(np.round(mean_ci_list[6]['zip'], 2)),''' & ''',np.round(mle_list[6]['beta-Poisson'][0], 2),tuple(np.round(ci_list[6]['beta-Poisson'][0], 2)),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['poisson'], 2),tuple(np.round(ci_list[7]['poisson'], 2)),''' & ''',np.round(mle_list[7]['geometric'], 2),tuple(np.round(ci_list[7]['geometric'], 2)),''' & ''',np.round(mle_list[7]['negative binomial'][0], 2),tuple(np.round(ci_list[7]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[7]['zip'][0]*(1-mle_list[7]['zip'][1]), 2),tuple(np.round(mean_ci_list[7]['zip'], 2)),''' & ''',np.round(mle_list[7]['beta-Poisson'][0], 2),tuple(np.round(ci_list[7]['beta-Poisson'][0], 2)),''' \\\\
    \\hline
	\end{tabular}
	\caption{Maximum likelihood estimates and confidence intervals (in parentheses) for mean of offspring distribution by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Observed} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(var_list[0]['sample'], 2),''' & ''',np.round(var_list[0]['geometric'], 2),tuple(np.round(var_ci_list[0]['geometric'], 2)),''' & ''',np.round(var_list[0]['negative binomial'], 2),tuple(np.round(var_ci_list[0]['negative binomial'], 2)),''' & ''',np.round(var_list[0]['zip'], 2),tuple(np.round(var_ci_list[0]['zip'], 2)),''' & ''',np.round(var_list[0]['beta-Poisson'], 2),tuple(np.round(var_ci_list[0]['beta-Poisson'], 2)),''' \\\\
		Mpox & ''',np.round(var_list[1]['sample'], 2),''' & ''',np.round(var_list[1]['geometric'], 2),tuple(np.round(var_ci_list[1]['geometric'], 2)),''' & ''',np.round(var_list[1]['negative binomial'], 2),tuple(np.round(var_ci_list[1]['negative binomial'], 2)),''' & ''',np.round(var_list[1]['zip'], 2),tuple(np.round(var_ci_list[1]['zip'], 2)),''' & ''',np.round(var_list[1]['beta-Poisson'], 2),tuple(np.round(var_ci_list[1]['beta-Poisson'], 2)),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(var_list[2]['sample'], 2),''' & ''',np.round(var_list[2]['geometric'], 2),tuple(np.round(var_ci_list[2]['geometric'], 2)),''' & ''',np.round(var_list[2]['negative binomial'], 2),tuple(np.round(var_ci_list[2]['negative binomial'], 2)),''' & ''',np.round(var_list[2]['zip'], 2),tuple(np.round(var_ci_list[2]['zip'], 2)),''' & ''',np.round(var_list[2]['beta-Poisson'], 2),tuple(np.round(var_ci_list[2]['beta-Poisson'], 2)),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(var_list[3]['sample'], 2),''' & ''',np.round(var_list[3]['geometric'], 2),tuple(np.round(var_ci_list[3]['geometric'], 2)),''' & ''',np.round(var_list[3]['negative binomial'], 2),tuple(np.round(var_ci_list[3]['negative binomial'], 2)),''' & ''',np.round(var_list[3]['zip'], 2),tuple(np.round(var_ci_list[3]['zip'], 2)),''' & ''',np.round(var_list[3]['beta-Poisson'], 2),tuple(np.round(var_ci_list[3]['beta-Poisson'], 2)),''' \\\\
		SARS, Singapore 2003 & ''',np.round(var_list[4]['sample'], 2),''' & ''',np.round(var_list[4]['geometric'], 2),tuple(np.round(var_ci_list[4]['geometric'], 2)),''' & ''',np.round(var_list[4]['negative binomial'], 2),tuple(np.round(var_ci_list[4]['negative binomial'], 2)),''' & ''',np.round(var_list[4]['zip'], 2),tuple(np.round(var_ci_list[4]['zip'], 2)),''' & ''',np.round(var_list[4]['beta-Poisson'], 2),tuple(np.round(var_ci_list[4]['beta-Poisson'], 2)),''' \\\\
		MERS, South Korea 2015 & ''',np.round(var_list[5]['sample'], 2),''' & ''',np.round(var_list[5]['geometric'], 2),tuple(np.round(var_ci_list[5]['geometric'], 2)),''' & ''',np.round(var_list[5]['negative binomial'], 2),tuple(np.round(var_ci_list[5]['negative binomial'], 2)),''' & ''',np.round(var_list[5]['zip'], 2),tuple(np.round(var_ci_list[5]['zip'], 2)),''' & ''',np.round(var_list[5]['beta-Poisson'], 2),tuple(np.round(var_ci_list[5]['beta-Poisson'], 2)),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(var_list[6]['sample'], 2),''' & ''',np.round(var_list[6]['geometric'], 2),tuple(np.round(var_ci_list[6]['geometric'], 2)),''' & ''',np.round(var_list[6]['negative binomial'], 2),tuple(np.round(var_ci_list[6]['negative binomial'], 2)),''' & ''',np.round(var_list[6]['zip'], 2),tuple(np.round(var_ci_list[6]['zip'], 2)),''' & ''',np.round(var_list[6]['beta-Poisson'], 2),tuple(np.round(var_ci_list[6]['beta-Poisson'], 2)),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(var_list[7]['sample'], 2),''' & ''',np.round(var_list[7]['geometric'], 2),tuple(np.round(var_ci_list[7]['geometric'], 2)),''' & ''',np.round(var_list[7]['negative binomial'], 2),tuple(np.round(var_ci_list[7]['negative binomial'], 2)),''' & ''',np.round(var_list[7]['zip'], 2),tuple(np.round(var_ci_list[7]['zip'], 2)),''' & ''',np.round(var_list[7]['beta-Poisson'], 2),tuple(np.round(var_ci_list[7]['beta-Poisson'], 2)),''' \\\\
    \\hline
	\end{tabular}
	\caption{Maximum likelihood estimates for variance of each maximum likelihood distribution, 95\\% confidence intervals in parentheses. The Poisson distribution is omitted since the MLE's and confidence intervals are identical to those of the mean.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Sample mean} $\\lambda$ & \\textbf{Neg. bin. } $\\theta$ & \\textbf{ZIP} $\\tilde{\lambda}$ & \\textbf{ZIP} $\sigma$ & \\textbf{Beta-Poisson} $\Phi$ & \\textbf{Beta-Poisson} $N$ \\
    \hline
		Plague & ''',np.round(mle_list[0]['poisson'], 2),tuple(np.round(ci_list[0]['poisson'], 2)),''' & ''',np.round(mle_list[0]['negative binomial'][1], 2),tuple(np.round(ci_list[0]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[0]['zip'][0], 2),tuple(np.round(ci_list[1]['zip'][0], 2)),''' & ''',np.round(mle_list[0]['zip'][1], 2),tuple(np.round(ci_list[0]['zip'][1], 2)),''' & ''',np.round(mle_list[0]['beta-Poisson'][1], 2),tuple(np.round(ci_list[0]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[0]['beta-Poisson'][2], 2),tuple(np.round(ci_list[0]['beta-Poisson'][2], 2)),'''\\
		Mpox & ''',np.round(mle_list[1]['poisson'], 2),tuple(np.round(ci_list[1]['poisson'], 2)),''' & ''',np.round(mle_list[1]['negative binomial'][1], 2),tuple(np.round(ci_list[1]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[1]['zip'][0], 2),tuple(np.round(ci_list[1]['zip'][0], 2)),''' & ''',np.round(mle_list[1]['zip'][1], 2),tuple(np.round(ci_list[1]['zip'][1], 2)),''' & ''',np.round(mle_list[1]['beta-Poisson'][1], 2),tuple(np.round(ci_list[1]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[1]['beta-Poisson'][2], 2),tuple(np.round(ci_list[1]['beta-Poisson'][2], 2)),'''\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['poisson'], 2),tuple(np.round(ci_list[2]['poisson'], 2)),''' & ''',np.round(mle_list[2]['negative binomial'][1], 2),tuple(np.round(ci_list[2]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[2]['zip'][0], 2),tuple(np.round(ci_list[2]['zip'][0], 2)),''' & ''',np.round(mle_list[2]['zip'][1], 2),tuple(np.round(ci_list[2]['zip'][1], 2)),''' & ''',np.round(mle_list[2]['beta-Poisson'][1], 2),tuple(np.round(ci_list[2]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[2]['beta-Poisson'][2], 2),tuple(np.round(ci_list[2]['beta-Poisson'][2], 2)),'''\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['poisson'], 2),tuple(np.round(ci_list[3]['poisson'], 2)),''' & ''',np.round(mle_list[3]['negative binomial'][1], 2),tuple(np.round(ci_list[3]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[3]['zip'][0], 2),tuple(np.round(ci_list[3]['zip'][0], 2)),''' & ''',np.round(mle_list[3]['zip'][1], 2),tuple(np.round(ci_list[3]['zip'][1], 2)),''' & ''',np.round(mle_list[3]['beta-Poisson'][1], 2),tuple(np.round(ci_list[3]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[3]['beta-Poisson'][2], 2),tuple(np.round(ci_list[3]['beta-Poisson'][2], 2)),'''\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['poisson'], 2),tuple(np.round(ci_list[4]['poisson'], 2)),''' & ''',np.round(mle_list[4]['negative binomial'][1], 2),tuple(np.round(ci_list[4]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[4]['zip'][0], 2),tuple(np.round(ci_list[4]['zip'][0], 2)),''' & ''',np.round(mle_list[4]['zip'][1], 2),tuple(np.round(ci_list[4]['zip'][1], 2)),''' & ''',np.round(mle_list[4]['beta-Poisson'][1], 2),tuple(np.round(ci_list[4]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[4]['beta-Poisson'][2], 2),tuple(np.round(ci_list[4]['beta-Poisson'][2], 2)),'''\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['poisson'], 2),tuple(np.round(ci_list[5]['poisson'], 2)),''' & ''',np.round(mle_list[5]['negative binomial'][1], 2),tuple(np.round(ci_list[5]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[5]['zip'][0], 2),tuple(np.round(ci_list[5]['zip'][0], 2)),''' & ''',np.round(mle_list[5]['zip'][1], 2),tuple(np.round(ci_list[5]['zip'][1], 2)),''' & ''',np.round(mle_list[5]['beta-Poisson'][1], 2),tuple(np.round(ci_list[5]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[5]['beta-Poisson'][2], 2),tuple(np.round(ci_list[5]['beta-Poisson'][2], 2)),'''\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['poisson'], 2),tuple(np.round(ci_list[6]['poisson'], 2)),''' & ''',np.round(mle_list[6]['negative binomial'][1], 2),tuple(np.round(ci_list[6]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[6]['zip'][0], 2),tuple(np.round(ci_list[6]['zip'][0], 2)),''' & ''',np.round(mle_list[6]['zip'][1], 2),tuple(np.round(ci_list[6]['zip'][1], 2)),''' & ''',np.round(mle_list[6]['beta-Poisson'][1], 2),tuple(np.round(ci_list[6]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[6]['beta-Poisson'][2], 2),tuple(np.round(ci_list[6]['beta-Poisson'][2], 2)),'''\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['poisson'], 2),tuple(np.round(ci_list[7]['poisson'], 2)),''' & ''',np.round(mle_list[7]['negative binomial'][1], 2),tuple(np.round(ci_list[7]['negative binomial'][1], 2)),''' & ''',np.round(mle_list[7]['zip'][0], 2),tuple(np.round(ci_list[7]['zip'][0], 2)),''' & ''',np.round(mle_list[7]['zip'][1], 2),tuple(np.round(ci_list[7]['zip'][1], 2)),''' & ''',np.round(mle_list[7]['beta-Poisson'][1], 2),tuple(np.round(ci_list[7]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[7]['beta-Poisson'][2], 2),tuple(np.round(ci_list[7]['beta-Poisson'][2], 2)),'''\\
	\\hline
	\end{tabular}
	\caption{Maximum likelihood estimates and confidence intervals (in parentheses) for model parameters given each dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & $\\lambda$ & $\\Phi$ & $\\nu$ \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['beta-Poisson'][0], 2),tuple(np.round(ci_list[0]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[0]['beta-Poisson'][1], 2),tuple(np.round(ci_list[0]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[0]['beta-Poisson'][2], 2),tuple(np.round(ci_list[0]['beta-Poisson'][2], 2)),'''\\\\
		Mpox & ''',np.round(mle_list[1]['beta-Poisson'][0], 2),tuple(np.round(ci_list[1]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[1]['beta-Poisson'][1], 2),tuple(np.round(ci_list[1]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[1]['beta-Poisson'][2], 2),tuple(np.round(ci_list[1]['beta-Poisson'][2], 2)),'''\\\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['beta-Poisson'][0], 2),tuple(np.round(ci_list[2]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[2]['beta-Poisson'][1], 2),tuple(np.round(ci_list[2]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[2]['beta-Poisson'][2], 2),tuple(np.round(ci_list[2]['beta-Poisson'][2], 2)),'''\\\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['beta-Poisson'][0], 2),tuple(np.round(ci_list[3]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[3]['beta-Poisson'][1], 2),tuple(np.round(ci_list[3]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[3]['beta-Poisson'][2], 2),tuple(np.round(ci_list[3]['beta-Poisson'][2], 2)),'''\\\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['beta-Poisson'][0], 2),tuple(np.round(ci_list[4]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[4]['beta-Poisson'][1], 2),tuple(np.round(ci_list[4]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[4]['beta-Poisson'][2], 2),tuple(np.round(ci_list[4]['beta-Poisson'][2], 2)),'''\\\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['beta-Poisson'][0], 2),tuple(np.round(ci_list[5]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[5]['beta-Poisson'][1], 2),tuple(np.round(ci_list[5]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[5]['beta-Poisson'][2], 2),tuple(np.round(ci_list[5]['beta-Poisson'][2], 2)),'''\\\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['beta-Poisson'][0], 2),tuple(np.round(ci_list[6]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[6]['beta-Poisson'][1], 2),tuple(np.round(ci_list[6]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[6]['beta-Poisson'][2], 2),tuple(np.round(ci_list[6]['beta-Poisson'][2], 2)),'''\\\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['beta-Poisson'][0], 2),tuple(np.round(ci_list[7]['beta-Poisson'][0], 2)),''' & ''',np.round(mle_list[7]['beta-Poisson'][1], 2),tuple(np.round(ci_list[7]['beta-Poisson'][1], 2)),''' & ''',np.round(mle_list[7]['beta-Poisson'][2], 2),tuple(np.round(ci_list[7]['beta-Poisson'][2], 2)),'''\\\\
	\\hline
	\end{tabular}
	\caption{Maximum likelihood estimates of beta-Poisson model parameters by dataset, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & $\\lambda$ $ \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['poisson'], 2),tuple(np.round(ci_list[0]['poisson'], 2)),'''\\\\
		Mpox & ''',np.round(mle_list[1]['poisson'], 2),tuple(np.round(ci_list[1]['poisson'], 2)),'''\\\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['poisson'], 2),tuple(np.round(ci_list[2]['poisson'], 2)),'''\\\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['poisson'], 2),tuple(np.round(ci_list[3]['poisson'], 2)),'''\\\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['poisson'], 2),tuple(np.round(ci_list[4]['poisson'], 2)),'''\\\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['poisson'], 2),tuple(np.round(ci_list[5]['poisson'], 2)),'''\\\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['poisson'], 2),tuple(np.round(ci_list[6]['poisson'], 2)),'''\\\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['poisson'], 2),tuple(np.round(ci_list[7]['poisson'], 2)),'''\\\\
	\\hline
	\end{tabular}
	\caption{Maximum likelihood estimates of Poisson model parameters by dataset, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & $\\lambda$ $ \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['geometric'], 2),tuple(np.round(ci_list[0]['geometric'], 2)),'''\\\\
		Mpox & ''',np.round(mle_list[1]['geometric'], 2),tuple(np.round(ci_list[1]['geometric'], 2)),'''\\\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['geometric'], 2),tuple(np.round(ci_list[2]['geometric'], 2)),'''\\\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['geometric'], 2),tuple(np.round(ci_list[3]['geometric'], 2)),'''\\\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['geometric'], 2),tuple(np.round(ci_list[4]['geometric'], 2)),'''\\\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['geometric'], 2),tuple(np.round(ci_list[5]['geometric'], 2)),'''\\\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['geometric'], 2),tuple(np.round(ci_list[6]['geometric'], 2)),'''\\\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['geometric'], 2),tuple(np.round(ci_list[7]['geometric'], 2)),'''\\\\
	\\hline
	\end{tabular}
	\caption{Maximum likelihood estimates of geometric model parameters by dataset, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & $\\lambda$ & $\\theta$ \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['negative binomial'][0], 2),tuple(np.round(ci_list[0]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[0]['negative binomial'][1], 2),tuple(np.round(ci_list[0]['negative binomial'][1], 2)),'''\\\\
        Mpox & ''',np.round(mle_list[1]['negative binomial'][0], 2),tuple(np.round(ci_list[1]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[1]['negative binomial'][1], 2),tuple(np.round(ci_list[1]['negative binomial'][1], 2)),'''\\\\
        Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['negative binomial'][0], 2),tuple(np.round(ci_list[2]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[2]['negative binomial'][1], 2),tuple(np.round(ci_list[2]['negative binomial'][1], 2)),'''\\\\
        Ebola, Guinea 2014 & ''',np.round(mle_list[3]['negative binomial'][0], 2),tuple(np.round(ci_list[3]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[3]['negative binomial'][1], 2),tuple(np.round(ci_list[3]['negative binomial'][1], 2)),'''\\\\
        SARS, Singapore 2003 & ''',np.round(mle_list[4]['negative binomial'][0], 2),tuple(np.round(ci_list[4]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[4]['negative binomial'][1], 2),tuple(np.round(ci_list[4]['negative binomial'][1], 2)),'''\\\\
        MERS, South Korea 2015 & ''',np.round(mle_list[5]['negative binomial'][0], 2),tuple(np.round(ci_list[5]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[5]['negative binomial'][1], 2),tuple(np.round(ci_list[5]['negative binomial'][1], 2)),'''\\\\
        MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['negative binomial'][0], 2),tuple(np.round(ci_list[6]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[6]['negative binomial'][1], 2),tuple(np.round(ci_list[6]['negative binomial'][1], 2)),'''\\\\
        Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['negative binomial'][0], 2),tuple(np.round(ci_list[7]['negative binomial'][0], 2)),''' & ''',np.round(mle_list[7]['negative binomial'][1], 2),tuple(np.round(ci_list[7]['negative binomial'][1], 2)),'''\\\\
    \hline
	\end{tabular}
	\caption{Maximum likelihood estimates of negative binomial model parameters by dataset, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & $\\tilde{\\lambda}$ & $\\sigma$ \\\\
    \hline
		Plague & ''',np.round(mle_list[0]['zip'][0], 2),tuple(np.round(ci_list[0]['zip'][0], 2)),''' & ''',np.round(mle_list[0]['zip'][1], 2),tuple(np.round(ci_list[0]['zip'][1], 2)),'''\\\\
        Mpox & ''',np.round(mle_list[1]['zip'][0], 2),tuple(np.round(ci_list[1]['zip'][0], 2)),''' & ''',np.round(mle_list[1]['zip'][1], 2),tuple(np.round(ci_list[1]['zip'][1], 2)),'''\\\\
        Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['zip'][0], 2),tuple(np.round(ci_list[2]['zip'][0], 2)),''' & ''',np.round(mle_list[2]['zip'][1], 2),tuple(np.round(ci_list[2]['zip'][1], 2)),'''\\\\
        Ebola, Guinea 2014 & ''',np.round(mle_list[3]['zip'][0], 2),tuple(np.round(ci_list[3]['zip'][0], 2)),''' & ''',np.round(mle_list[3]['zip'][1], 2),tuple(np.round(ci_list[3]['zip'][1], 2)),'''\\\\
        SARS, Singapore 2003 & ''',np.round(mle_list[4]['zip'][0], 2),tuple(np.round(ci_list[4]['zip'][0], 2)),''' & ''',np.round(mle_list[4]['zip'][1], 2),tuple(np.round(ci_list[4]['zip'][1], 2)),'''\\\\
        MERS, South Korea 2015 & ''',np.round(mle_list[5]['zip'][0], 2),tuple(np.round(ci_list[5]['zip'][0], 2)),''' & ''',np.round(mle_list[5]['zip'][1], 2),tuple(np.round(ci_list[5]['zip'][1], 2)),'''\\\\
        MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['zip'][0], 2),tuple(np.round(ci_list[6]['zip'][0], 2)),''' & ''',np.round(mle_list[6]['zip'][1], 2),tuple(np.round(ci_list[6]['zip'][1], 2)),'''\\\\
        Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['zip'][0], 2),tuple(np.round(ci_list[7]['zip'][0], 2)),''' & ''',np.round(mle_list[7]['zip'][1], 2),tuple(np.round(ci_list[7]['zip'][1], 2)),'''\\\\
    \hline
	\end{tabular}
	\caption{Maximum likelihood estimates of zero-inflated Poisson model parameters by dataset, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lc}
		\\hline
		\\textbf{Dataset} & $\\bar{\\nu}$ \\
    \hline
		Plague & ''',np.round(mle_list[0]['beta-Poisson'][2], 2),tuple(np.round(ci_list[0]['beta-Poisson'][2], 2)),'''\\\\
		Mpox & ''',np.round(mle_list[1]['beta-Poisson'][2], 2),tuple(np.round(ci_list[1]['beta-Poisson'][2], 2)),'''\\\\
		Ebola, Nigeria 2014 & ''',np.round(mle_list[2]['beta-Poisson'][2], 2),tuple(np.round(ci_list[2]['beta-Poisson'][2], 2)),'''\\\\
		Ebola, Guinea 2014 & ''',np.round(mle_list[3]['beta-Poisson'][2], 2),tuple(np.round(ci_list[3]['beta-Poisson'][2], 2)),'''\\\\
		SARS, Singapore 2003 & ''',np.round(mle_list[4]['beta-Poisson'][2], 2),tuple(np.round(ci_list[4]['beta-Poisson'][2], 2)),'''\\\\
		MERS, South Korea 2015 & ''',np.round(mle_list[5]['beta-Poisson'][2], 2),tuple(np.round(ci_list[5]['beta-Poisson'][2], 2)),'''\\\\
		MERS, Saudi Arabia 2015 & ''',np.round(mle_list[6]['beta-Poisson'][2], 2),tuple(np.round(ci_list[6]['beta-Poisson'][2], 2)),'''\\\\
		Norovirus, Netherlands 2012 & ''',np.round(mle_list[7]['beta-Poisson'][2], 2),tuple(np.round(ci_list[7]['beta-Poisson'][2], 2)),'''\\\\
	\\hline
	\end{tabular}
	\caption{Maximum likelihood estimates of inverse contact parameter $\\nu$ for each dataset.}
\end{table}
''')



print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(llh_list[0]['poisson'][0], 2),''' & ''',np.round(llh_list[0]['geometric'][0], 2),''' & ''',np.round(llh_list[0]['negative binomial'][0], 2),''' & ''',np.round(llh_list[0]['zip'][0], 2),''' & ''',np.round(llh_list[0]['beta-Poisson'][0], 2),''' \\\\
		Mpox & ''',np.round(llh_list[1]['poisson'][0], 2),''' & ''',np.round(llh_list[1]['geometric'][0], 2),''' & ''',np.round(llh_list[1]['negative binomial'][0], 2),''' & ''',np.round(llh_list[1]['zip'][0], 2),''' & ''',np.round(llh_list[1]['beta-Poisson'][0], 2),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(llh_list[2]['poisson'][0], 2),''' & ''',np.round(llh_list[2]['geometric'][0], 2),''' & ''',np.round(llh_list[2]['negative binomial'][0], 2),''' & ''',np.round(llh_list[2]['zip'][0], 2),''' & ''',np.round(llh_list[2]['beta-Poisson'][0], 2),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(llh_list[3]['poisson'][0], 2),''' & ''',np.round(llh_list[3]['geometric'][0], 2),''' & ''',np.round(llh_list[3]['negative binomial'][0], 2),''' & ''',np.round(llh_list[3]['zip'][0], 2),''' & ''',np.round(llh_list[3]['beta-Poisson'][0], 2),''' \\\\
		SARS, Singapore 2003 & ''',np.round(llh_list[4]['poisson'][0], 2),''' & ''',np.round(llh_list[4]['geometric'][0], 2),''' & ''',np.round(llh_list[4]['negative binomial'][0], 2),''' & ''',np.round(llh_list[4]['zip'][0], 2),''' & ''',np.round(llh_list[4]['beta-Poisson'][0], 2),''' \\\\
		MERS, South Korea 2015 & ''',np.round(llh_list[5]['poisson'][0], 2),''' & ''',np.round(llh_list[5]['geometric'][0], 2),''' & ''',np.round(llh_list[5]['negative binomial'][0], 2),''' & ''',np.round(llh_list[5]['zip'][0], 2),''' & ''',np.round(llh_list[5]['beta-Poisson'][0], 2),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(llh_list[6]['poisson'][0], 2),''' & ''',np.round(llh_list[6]['geometric'][0], 2),''' & ''',np.round(llh_list[6]['negative binomial'][0], 2),''' & ''',np.round(llh_list[6]['zip'][0], 2),''' & ''',np.round(llh_list[6]['beta-Poisson'][0], 2),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(llh_list[7]['poisson'][0], 2),''' & ''',np.round(llh_list[7]['geometric'][0], 2),''' & ''',np.round(llh_list[7]['negative binomial'][0], 2),''' & ''',np.round(llh_list[7]['zip'][0], 2),''' & ''',np.round(llh_list[7]['beta-Poisson'][0], 2),''' \\\\
    \\hline
	\end{tabular}
	\caption{Log likelihood at MLE by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(llh_list[0]['poisson'][1], 2),''' & ''',np.round(llh_list[0]['geometric'][1], 2),''' & ''',np.round(llh_list[0]['negative binomial'][1], 2),''' & ''',np.round(llh_list[0]['zip'][1], 2),''' & ''',np.round(llh_list[0]['beta-Poisson'][1], 2),''' \\\\
		Mpox & ''',np.round(llh_list[1]['poisson'][1], 2),''' & ''',np.round(llh_list[1]['geometric'][1], 2),''' & ''',np.round(llh_list[1]['negative binomial'][1], 2),''' & ''',np.round(llh_list[1]['zip'][1], 2),''' & ''',np.round(llh_list[1]['beta-Poisson'][1], 2),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(llh_list[2]['poisson'][1], 2),''' & ''',np.round(llh_list[2]['geometric'][1], 2),''' & ''',np.round(llh_list[2]['negative binomial'][1], 2),''' & ''',np.round(llh_list[2]['zip'][1], 2),''' & ''',np.round(llh_list[2]['beta-Poisson'][1], 2),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(llh_list[3]['poisson'][1], 2),''' & ''',np.round(llh_list[3]['geometric'][1], 2),''' & ''',np.round(llh_list[3]['negative binomial'][1], 2),''' & ''',np.round(llh_list[3]['zip'][1], 2),''' & ''',np.round(llh_list[3]['beta-Poisson'][1], 2),''' \\\\
		SARS, Singapore 2003 & ''',np.round(llh_list[4]['poisson'][1], 2),''' & ''',np.round(llh_list[4]['geometric'][1], 2),''' & ''',np.round(llh_list[4]['negative binomial'][1], 2),''' & ''',np.round(llh_list[4]['zip'][1], 2),''' & ''',np.round(llh_list[4]['beta-Poisson'][1], 2),''' \\\\
		MERS, South Korea 2015 & ''',np.round(llh_list[5]['poisson'][1], 2),''' & ''',np.round(llh_list[5]['geometric'][1], 2),''' & ''',np.round(llh_list[5]['negative binomial'][1], 2),''' & ''',np.round(llh_list[5]['zip'][1], 2),''' & ''',np.round(llh_list[5]['beta-Poisson'][1], 2),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(llh_list[6]['poisson'][1], 2),''' & ''',np.round(llh_list[6]['geometric'][1], 2),''' & ''',np.round(llh_list[6]['negative binomial'][1], 2),''' & ''',np.round(llh_list[6]['zip'][1], 2),''' & ''',np.round(llh_list[6]['beta-Poisson'][1], 2),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(llh_list[7]['poisson'][1], 2),''' & ''',np.round(llh_list[7]['geometric'][1], 2),''' & ''',np.round(llh_list[7]['negative binomial'][1], 2),''' & ''',np.round(llh_list[7]['zip'][1], 2),''' & ''',np.round(llh_list[7]['beta-Poisson'][1], 2),''' \\\\
    \\hline
	\end{tabular}
	\caption{Akaike information criterion at MLE by model and dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} \\\\
    \hline
		Plague & ''',np.round(np.exp(llh_list[0]['poisson'][0]-llh_list[0]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[0]['geometric'][0]-llh_list[0]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[0]['negative binomial'][0]-llh_list[0]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[0]['zip'][0]-llh_list[0]['beta-Poisson'][0]), 2),''' \\\\
		Mpox & ''',np.round(np.exp(llh_list[1]['poisson'][0]-llh_list[1]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[1]['geometric'][0]-llh_list[1]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[1]['negative binomial'][0]-llh_list[1]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[1]['zip'][0]-llh_list[1]['beta-Poisson'][0]), 2),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(np.exp(llh_list[2]['poisson'][0]-llh_list[2]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[2]['geometric'][0]-llh_list[2]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[2]['negative binomial'][0]-llh_list[2]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[2]['zip'][0]-llh_list[2]['beta-Poisson'][0]), 2),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(np.exp(llh_list[3]['poisson'][0]-llh_list[3]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[3]['geometric'][0]-llh_list[3]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[3]['negative binomial'][0]-llh_list[3]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[3]['zip'][0]-llh_list[3]['beta-Poisson'][0]), 2),''' \\\\
		SARS, Singapore 2003 & ''',np.round(np.exp(llh_list[4]['poisson'][0]-llh_list[4]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[4]['geometric'][0]-llh_list[4]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[4]['negative binomial'][0]-llh_list[4]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[4]['zip'][0]-llh_list[4]['beta-Poisson'][0]), 2),''' \\\\
		MERS, South Korea 2015 & ''',np.round(np.exp(llh_list[5]['poisson'][0]-llh_list[5]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[5]['geometric'][0]-llh_list[5]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[5]['negative binomial'][0]-llh_list[5]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[5]['zip'][0]-llh_list[5]['beta-Poisson'][0]), 2),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(np.exp(llh_list[6]['poisson'][0]-llh_list[6]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[6]['geometric'][0]-llh_list[6]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[6]['negative binomial'][0]-llh_list[6]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[6]['zip'][0]-llh_list[6]['beta-Poisson'][0]), 2),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(np.exp(llh_list[7]['poisson'][0]-llh_list[7]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[7]['geometric'][0]-llh_list[7]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[7]['negative binomial'][0]-llh_list[7]['beta-Poisson'][0]), 2),''' & ''',np.round(np.exp(llh_list[7]['zip'][0]-llh_list[7]['beta-Poisson'][0]), 2),''' \\\\
	\\hline
	\end{tabular}
	\caption{Likelihood ratios of beta-Poisson to negative binomial and ZIP models under each dataset.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} \\\\
    \hline
		Plague & ''',np.round(chi2.sf(-2*(llh_list[0]['poisson'][0]-llh_list[0]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[0]['geometric'][0]-llh_list[0]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[0]['negative binomial'][0]-llh_list[0]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[0]['zip'][0]-llh_list[0]['beta-Poisson'][0]), 1), 2),''' \\\\
		Mpox & ''',np.round(chi2.sf(-2*(llh_list[1]['poisson'][0]-llh_list[1]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[1]['geometric'][0]-llh_list[1]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[1]['negative binomial'][0]-llh_list[1]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[1]['zip'][0]-llh_list[1]['beta-Poisson'][0]), 1), 2),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(chi2.sf(-2*(llh_list[2]['poisson'][0]-llh_list[2]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[2]['geometric'][0]-llh_list[2]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[2]['negative binomial'][0]-llh_list[2]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[2]['zip'][0]-llh_list[2]['beta-Poisson'][0]), 1), 2),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(chi2.sf(-2*(llh_list[3]['poisson'][0]-llh_list[3]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[3]['geometric'][0]-llh_list[3]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[3]['negative binomial'][0]-llh_list[3]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[3]['zip'][0]-llh_list[3]['beta-Poisson'][0]), 1), 2),''' \\\\
		SARS, Singapore 2003 & ''',np.round(chi2.sf(-2*(llh_list[4]['poisson'][0]-llh_list[4]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[4]['geometric'][0]-llh_list[4]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[4]['negative binomial'][0]-llh_list[4]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[4]['zip'][0]-llh_list[4]['beta-Poisson'][0]), 1), 2),''' \\\\
		MERS, South Korea 2015 & ''',np.round(chi2.sf(-2*(llh_list[5]['poisson'][0]-llh_list[5]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[5]['geometric'][0]-llh_list[5]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[5]['negative binomial'][0]-llh_list[5]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[5]['zip'][0]-llh_list[5]['beta-Poisson'][0]), 1), 2),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(chi2.sf(-2*(llh_list[6]['poisson'][0]-llh_list[6]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[6]['geometric'][0]-llh_list[6]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[6]['negative binomial'][0]-llh_list[6]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[6]['zip'][0]-llh_list[6]['beta-Poisson'][0]), 1), 2),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(chi2.sf(-2*(llh_list[7]['poisson'][0]-llh_list[7]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[7]['geometric'][0]-llh_list[7]['beta-Poisson'][0]), 2), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[7]['negative binomial'][0]-llh_list[7]['beta-Poisson'][0]), 1), 2),''' & ''',np.round(chi2.sf(-2*(llh_list[7]['zip'][0]-llh_list[7]['beta-Poisson'][0]), 1), 2),''' \\\\
	\\hline
	\end{tabular}
	\caption{Likelihood ratio test statistics obtained by comparing beta-Poisson to other candidate models under each dataset.}
\end{table}
''')


print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Boundary} & \\textbf{Observed} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',superspread_list[0]['boundary'],''' & ''',np.round(superspread_list[0]['sample'], 2),''' & ''',np.round(superspread_list[0]['geometric'], 2),tuple(np.round(superspread_ci_list[0]['geometric'], 2)),''' & ''',np.round(superspread_list[0]['negative binomial'], 2),tuple(np.round(superspread_ci_list[0]['negative binomial'], 2)),''' & ''',np.round(superspread_list[0]['zip'], 2),tuple(np.round(superspread_ci_list[0]['zip'], 2)),''' & ''',np.round(superspread_list[0]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[0]['beta-Poisson'], 2)),''' \\\\
		Mpox & ''',superspread_list[1]['boundary'],''' & ''',np.round(superspread_list[1]['sample'], 2),''' & ''',np.round(superspread_list[1]['geometric'], 2),tuple(np.round(superspread_ci_list[1]['geometric'], 2)),''' & ''',np.round(superspread_list[1]['negative binomial'], 2),tuple(np.round(superspread_ci_list[1]['negative binomial'], 2)),''' & ''',np.round(superspread_list[1]['zip'], 2),tuple(np.round(superspread_ci_list[1]['zip'], 2)),''' & ''',np.round(superspread_list[1]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[1]['beta-Poisson'], 2)),''' \\\\
		Ebola, Nigeria 2014 & ''',superspread_list[2]['boundary'],''' & ''',np.round(superspread_list[2]['sample'], 2),''' & ''',np.round(superspread_list[2]['geometric'], 2),tuple(np.round(superspread_ci_list[2]['geometric'], 2)),''' & ''',np.round(superspread_list[2]['negative binomial'], 2),tuple(np.round(superspread_ci_list[2]['negative binomial'], 2)),''' & ''',np.round(superspread_list[2]['zip'], 2),tuple(np.round(superspread_ci_list[2]['zip'], 2)),''' & ''',np.round(superspread_list[2]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[2]['beta-Poisson'], 2)),''' \\\\
		Ebola, Guinea 2014 & ''',superspread_list[3]['boundary'],''' & ''',np.round(superspread_list[3]['sample'], 2),''' & ''',np.round(superspread_list[3]['geometric'], 2),tuple(np.round(superspread_ci_list[3]['geometric'], 2)),''' & ''',np.round(superspread_list[3]['negative binomial'], 2),tuple(np.round(superspread_ci_list[3]['negative binomial'], 2)),''' & ''',np.round(superspread_list[3]['zip'], 2),tuple(np.round(superspread_ci_list[3]['zip'], 2)),''' & ''',np.round(superspread_list[3]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[3]['beta-Poisson'], 2)),''' \\\\
		SARS, Singapore 2003 & ''',superspread_list[4]['boundary'],''' & ''',np.round(superspread_list[4]['sample'], 2),''' & ''',np.round(superspread_list[4]['geometric'], 2),tuple(np.round(superspread_ci_list[4]['geometric'], 2)),''' & ''',np.round(superspread_list[4]['negative binomial'], 2),tuple(np.round(superspread_ci_list[4]['negative binomial'], 2)),''' & ''',np.round(superspread_list[4]['zip'], 2),tuple(np.round(superspread_ci_list[4]['zip'], 2)),''' & ''',np.round(superspread_list[4]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[4]['beta-Poisson'], 2)),''' \\\\
		MERS, South Korea 2015 & ''',superspread_list[5]['boundary'],''' & ''',np.round(superspread_list[5]['sample'], 2),''' & ''',np.round(superspread_list[5]['geometric'], 2),tuple(np.round(superspread_ci_list[5]['geometric'], 2)),''' & ''',np.round(superspread_list[5]['negative binomial'], 2),tuple(np.round(superspread_ci_list[5]['negative binomial'], 2)),''' & ''',np.round(superspread_list[5]['zip'], 2),tuple(np.round(superspread_ci_list[5]['zip'], 2)),''' & ''',np.round(superspread_list[5]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[5]['beta-Poisson'], 2)),''' \\\\
		MERS, Saudi Arabia 2015 & ''',superspread_list[6]['boundary'],''' & ''',np.round(superspread_list[6]['sample'], 2),''' & ''',np.round(superspread_list[6]['geometric'], 2),tuple(np.round(superspread_ci_list[6]['geometric'], 2)),''' & ''',np.round(superspread_list[6]['negative binomial'], 2),tuple(np.round(superspread_ci_list[6]['negative binomial'], 2)),''' & ''',np.round(superspread_list[6]['zip'], 2),tuple(np.round(superspread_ci_list[6]['zip'], 2)),''' & ''',np.round(superspread_list[6]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[6]['beta-Poisson'], 2)),''' \\\\
		Norovirus, Netherlands 2012 & ''',superspread_list[7]['boundary'],''' & ''',np.round(superspread_list[7]['sample'], 2),''' & ''',np.round(superspread_list[7]['geometric'], 2),tuple(np.round(superspread_ci_list[7]['geometric'], 2)),''' & ''',np.round(superspread_list[7]['negative binomial'], 2),tuple(np.round(superspread_ci_list[7]['negative binomial'], 2)),''' & ''',np.round(superspread_list[7]['zip'], 2),tuple(np.round(superspread_ci_list[7]['zip'], 2)),''' & ''',np.round(superspread_list[7]['beta-Poisson'], 2),tuple(np.round(superspread_ci_list[7]['beta-Poisson'], 2)),''' \\\\
    \\hline
	\end{tabular}
	\caption{Superspreading boundary (99th percentile of fitted Poisson distribution) and proportion of cases above this boundary for each maximum likelihood distribution, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Poisson} & \\textbf{Observed} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(p0_list[0]['sample'], 2),''' & ''',np.round(p0_list[0]['poisson'], 2),tuple(np.round(p0_ci_list[0]['poisson'], 2)),''' & ''',np.round(p0_list[0]['geometric'], 2),tuple(np.round(p0_ci_list[0]['geometric'], 2)),''' & ''',np.round(p0_list[0]['negative binomial'], 2),tuple(np.round(p0_ci_list[0]['negative binomial'], 2)),''' & ''',np.round(p0_list[0]['zip'], 2),tuple(np.round(p0_ci_list[0]['zip'], 2)),''' & ''',np.round(p0_list[0]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[0]['beta-Poisson'], 2)),''' \\\\
		Mpox & ''',np.round(p0_list[1]['sample'], 2),''' & ''',np.round(p0_list[1]['poisson'], 2),tuple(np.round(p0_ci_list[1]['poisson'], 2)),''' & ''',np.round(p0_list[1]['geometric'], 2),tuple(np.round(p0_ci_list[1]['geometric'], 2)),''' & ''',np.round(p0_list[1]['negative binomial'], 2),tuple(np.round(p0_ci_list[1]['negative binomial'], 2)),''' & ''',np.round(p0_list[1]['zip'], 2),tuple(np.round(p0_ci_list[1]['zip'], 2)),''' & ''',np.round(p0_list[1]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[1]['beta-Poisson'], 2)),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(p0_list[2]['sample'], 2),''' & ''',np.round(p0_list[2]['poisson'], 2),tuple(np.round(p0_ci_list[2]['poisson'], 2)),''' & ''',np.round(p0_list[2]['geometric'], 2),tuple(np.round(p0_ci_list[2]['geometric'], 2)),''' & ''',np.round(p0_list[2]['negative binomial'], 2),tuple(np.round(p0_ci_list[2]['negative binomial'], 2)),''' & ''',np.round(p0_list[2]['zip'], 2),tuple(np.round(p0_ci_list[2]['zip'], 2)),''' & ''',np.round(p0_list[2]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[2]['beta-Poisson'], 2)),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(p0_list[3]['sample'], 2),''' & ''',np.round(p0_list[3]['poisson'], 2),tuple(np.round(p0_ci_list[3]['poisson'], 2)),''' & ''',np.round(p0_list[3]['geometric'], 2),tuple(np.round(p0_ci_list[3]['geometric'], 2)),''' & ''',np.round(p0_list[3]['negative binomial'], 2),tuple(np.round(p0_ci_list[3]['negative binomial'], 2)),''' & ''',np.round(p0_list[3]['zip'], 2),tuple(np.round(p0_ci_list[3]['zip'], 2)),''' & ''',np.round(p0_list[3]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[3]['beta-Poisson'], 2)),''' \\\\
		SARS, Singapore 2003 & ''',np.round(p0_list[4]['sample'], 2),''' & ''',np.round(p0_list[4]['poisson'], 2),tuple(np.round(p0_ci_list[4]['poisson'], 2)),''' & ''',np.round(p0_list[4]['geometric'], 2),tuple(np.round(p0_ci_list[4]['geometric'], 2)),''' & ''',np.round(p0_list[4]['negative binomial'], 2),tuple(np.round(p0_ci_list[4]['negative binomial'], 2)),''' & ''',np.round(p0_list[4]['zip'], 2),tuple(np.round(p0_ci_list[4]['zip'], 2)),''' & ''',np.round(p0_list[4]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[4]['beta-Poisson'], 2)),''' \\\\
		MERS, South Korea 2015 & ''',np.round(p0_list[5]['sample'], 2),''' & ''',np.round(p0_list[5]['poisson'], 2),tuple(np.round(p0_ci_list[5]['poisson'], 2)),''' & ''',np.round(p0_list[5]['geometric'], 2),tuple(np.round(p0_ci_list[5]['geometric'], 2)),''' & ''',np.round(p0_list[5]['negative binomial'], 2),tuple(np.round(p0_ci_list[5]['negative binomial'], 2)),''' & ''',np.round(p0_list[5]['zip'], 2),tuple(np.round(p0_ci_list[5]['zip'], 2)),''' & ''',np.round(p0_list[5]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[5]['beta-Poisson'], 2)),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(p0_list[6]['sample'], 2),''' & ''',np.round(p0_list[6]['poisson'], 2),tuple(np.round(p0_ci_list[6]['poisson'], 2)),''' & ''',np.round(p0_list[6]['geometric'], 2),tuple(np.round(p0_ci_list[6]['geometric'], 2)),''' & ''',np.round(p0_list[6]['negative binomial'], 2),tuple(np.round(p0_ci_list[6]['negative binomial'], 2)),''' & ''',np.round(p0_list[6]['zip'], 2),tuple(np.round(p0_ci_list[6]['zip'], 2)),''' & ''',np.round(p0_list[6]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[6]['beta-Poisson'], 2)),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(p0_list[7]['sample'], 2),''' & ''',np.round(p0_list[7]['poisson'], 2),tuple(np.round(p0_ci_list[7]['poisson'], 2)),''' & ''',np.round(p0_list[7]['geometric'], 2),tuple(np.round(p0_ci_list[7]['geometric'], 2)),''' & ''',np.round(p0_list[7]['negative binomial'], 2),tuple(np.round(p0_ci_list[7]['negative binomial'], 2)),''' & ''',np.round(p0_list[7]['zip'], 2),tuple(np.round(p0_ci_list[7]['zip'], 2)),''' & ''',np.round(p0_list[7]['beta-Poisson'], 2),tuple(np.round(p0_ci_list[7]['beta-Poisson'], 2)),''' \\\\
    \\hline
	\end{tabular}
	\caption{Probabliity of a case generating zero secondary cases under each maximum likelihood distribution, 95\\% confidence intervals in parentheses.}
\end{table}
''')

print('''
\\begin{table}[ht]
	\centering
	\\begin{tabular}{lcccccc}
		\\hline
		\\textbf{Dataset} & \\textbf{Observed} & \\textbf{Geometric} & \\textbf{Neg. bin. } & \\textbf{ZIP} & \\textbf{Beta-Poisson} \\\\
    \hline
		Plague & ''',np.round(od_list[0]['sample'], 2),''' & ''',np.round(od_list[0]['geometric'], 2),tuple(np.round(od_ci_list[0]['geometric'], 2)),''' & ''',np.round(od_list[0]['negative binomial'], 2),tuple(np.round(od_ci_list[0]['negative binomial'], 2)),''' & ''',np.round(od_list[0]['zip'], 2),tuple(np.round(od_ci_list[0]['zip'], 2)),''' & ''',np.round(od_list[0]['beta-Poisson'], 2),tuple(np.round(od_ci_list[0]['beta-Poisson'], 2)),''' \\\\
		Mpox & ''',np.round(od_list[1]['sample'], 2),''' & ''',np.round(od_list[1]['geometric'], 2),tuple(np.round(od_ci_list[1]['geometric'], 2)),''' & ''',np.round(od_list[1]['negative binomial'], 2),tuple(np.round(od_ci_list[1]['negative binomial'], 2)),''' & ''',np.round(od_list[1]['zip'], 2),tuple(np.round(od_ci_list[1]['zip'], 2)),''' & ''',np.round(od_list[1]['beta-Poisson'], 2),tuple(np.round(od_ci_list[1]['beta-Poisson'], 2)),''' \\\\
		Ebola, Nigeria 2014 & ''',np.round(od_list[2]['sample'], 2),''' & ''',np.round(od_list[2]['geometric'], 2),tuple(np.round(od_ci_list[2]['geometric'], 2)),''' & ''',np.round(od_list[2]['negative binomial'], 2),tuple(np.round(od_ci_list[2]['negative binomial'], 2)),''' & ''',np.round(od_list[2]['zip'], 2),tuple(np.round(od_ci_list[2]['zip'], 2)),''' & ''',np.round(od_list[2]['beta-Poisson'], 2),tuple(np.round(od_ci_list[2]['beta-Poisson'], 2)),''' \\\\
		Ebola, Guinea 2014 & ''',np.round(od_list[3]['sample'], 2),''' & ''',np.round(od_list[3]['geometric'], 2),tuple(np.round(od_ci_list[3]['geometric'], 2)),''' & ''',np.round(od_list[3]['negative binomial'], 2),tuple(np.round(od_ci_list[3]['negative binomial'], 2)),''' & ''',np.round(od_list[3]['zip'], 2),tuple(np.round(od_ci_list[3]['zip'], 2)),''' & ''',np.round(od_list[3]['beta-Poisson'], 2),tuple(np.round(od_ci_list[3]['beta-Poisson'], 2)),''' \\\\
		SARS, Singapore 2003 & ''',np.round(od_list[4]['sample'], 2),''' & ''',np.round(od_list[4]['geometric'], 2),tuple(np.round(od_ci_list[4]['geometric'], 2)),''' & ''',np.round(od_list[4]['negative binomial'], 2),tuple(np.round(od_ci_list[4]['negative binomial'], 2)),''' & ''',np.round(od_list[4]['zip'], 2),tuple(np.round(od_ci_list[4]['zip'], 2)),''' & ''',np.round(od_list[4]['beta-Poisson'], 2),tuple(np.round(od_ci_list[4]['beta-Poisson'], 2)),''' \\\\
		MERS, South Korea 2015 & ''',np.round(od_list[5]['sample'], 2),''' & ''',np.round(od_list[5]['geometric'], 2),tuple(np.round(od_ci_list[5]['geometric'], 2)),''' & ''',np.round(od_list[5]['negative binomial'], 2),tuple(np.round(od_ci_list[5]['negative binomial'], 2)),''' & ''',np.round(od_list[5]['zip'], 2),tuple(np.round(od_ci_list[5]['zip'], 2)),''' & ''',np.round(od_list[5]['beta-Poisson'], 2),tuple(np.round(od_ci_list[5]['beta-Poisson'], 2)),''' \\\\
		MERS, Saudi Arabia 2015 & ''',np.round(od_list[6]['sample'], 2),''' & ''',np.round(od_list[6]['geometric'], 2),tuple(np.round(od_ci_list[6]['geometric'], 2)),''' & ''',np.round(od_list[6]['negative binomial'], 2),tuple(np.round(od_ci_list[6]['negative binomial'], 2)),''' & ''',np.round(od_list[6]['zip'], 2),tuple(np.round(od_ci_list[6]['zip'], 2)),''' & ''',np.round(od_list[6]['beta-Poisson'], 2),tuple(np.round(od_ci_list[6]['beta-Poisson'], 2)),''' \\\\
		Norovirus, Netherlands 2012 & ''',np.round(od_list[7]['sample'], 2),''' & ''',np.round(od_list[7]['geometric'], 2),tuple(np.round(od_ci_list[7]['geometric'], 2)),''' & ''',np.round(od_list[7]['negative binomial'], 2),tuple(np.round(od_ci_list[7]['negative binomial'], 2)),''' & ''',np.round(od_list[7]['zip'], 2),tuple(np.round(od_ci_list[7]['zip'], 2)),''' & ''',np.round(od_list[7]['beta-Poisson'], 2),tuple(np.round(od_ci_list[7]['beta-Poisson'], 2)),''' \\\\
    \\hline
	\end{tabular}
	\caption{Overdispersion of each maximum likelihood distribution, 95\\% confidence intervals in parentheses.}
\end{table}
''')