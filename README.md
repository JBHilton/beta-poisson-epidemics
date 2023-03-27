# beta-poisson-epidemic
The notebook(s) in this repository are designed to accompany an upcoming paper by Joe Hilton and Ian Hall which introduces a beta-Poisson mixture model for secondary case count data in the early stages of an infectious disease outbreak.

The secondary case datasets used in our study are stored in `datasets.py`.

The functions used in our study are stored in `functions.py`. These include functions for:
* calculating the log likelihood of parameters given data for each of the candidate models from our study (Poisson, geometric, negative binomial, zero-inflated Poisson (ZIP), and beta-Poisson);
* calculating maximum likelihood estimates for the parameters of the negative binomial, ZIP, and beta-Poisson distributions;
* calculating the probability mass function of the ZIP and beta-Poisson models;
* calculating the probability generating function of each of the candidate models;
* calculating the extinction probability of a branching process with offspring distribution given by each of the candidate models;
* creating a dictionary object given some data, which stores the parameter MLE's for each candidate model;
* creating a dictionary object which stores the MLE's of the variance of each candidate model given some data;
* creating a dictionary object which stores the MLE's of the superspreading proportion of each candidate model given some data;
* creating a dictionary object which stores the MLE's of the proportion of zeros in each candidate model given some data;
* creating a dictionary object which stores the log likelihoods and AIC's of each candidate model with maximum likelihood parameters given some data;
* calculates confidence intervals for a statistic from a set of bootstrap samples.

The script `parallel_confidence_intervals.py` calculates MLEs and confidence intervals of model parameters and statistics for each candidate model for a dataset chosen by the user from those stored in `datasets.py`, with the number of bootstraps used in the confidence interval calculation also chosen by the user. The script `sensitivity_analyses.py` performs a sensitivity analysis for a dataset chosen by the user, with log likelihoods given the dataset calculated over a range of beta-Poisson model parameters. The script `plot_outputs.py` generates the plots presented in our manuscript from the outputs of `parallel_confidence_intervals.py`, while `plot_llh_curves.py` generates the plot presented in our appendix from the outputs of `sensititivity_analysis.py`.
