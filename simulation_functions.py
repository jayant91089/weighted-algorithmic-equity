import numpy as np
import pandas as pd

def expit(x):
    return np.exp(x) / (1 + np.exp(x))

def p_trunc(probs):
    return np.clip(probs, 0.005, 0.995)

# Variables returned by calc_cov_bias
return_vars_covbias = [
    "fnr_true", "fnr_wgt", "tpr_true", "tpr_wgt", "fpr_true", "fpr_wgt", "tnr_true", "tnr_wgt", "ppv_true", "ppv_wgt", "npv_true", "npv_wgt",
    "bias_obs_fnr", "bias_obs_abs_fnr", "bias_obs_tpr", "bias_obs_abs_tpr", "bias_obs_fpr", "bias_obs_abs_fpr", "bias_obs_tnr", "bias_obs_abs_tnr", 
    "bias_obs_ppv", "bias_obs_abs_ppv", "bias_obs_npv", "bias_obs_abs_npv",
    "fnr_epsilon", "tpr_epsilon", "fpr_epsilon", "tnr_epsilon", "ppv_epsilon", "npv_epsilon",
    "pa_ratio", "bias_true_fnr",
    "bias_bound_abs_fnr", "bias_bound_abs_tpr", "bias_bound_abs_fpr", "bias_bound_abs_tnr", "bias_bound_abs_ppv", "bias_bound_abs_npv",
    "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr", "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
    "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr", "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv",
    "bias_epsilon_fnr", "bias_epsilon_tpr", "bias_epsilon_fpr", "bias_epsilon_tnr", "bias_epsilon_ppv", "bias_epsilon_npv"
]

# Helper function to get bias terms across multiple A groups
def get_avals_bias(dat, epsilon, epsilon_prime):
    avals = dat['A'].unique()
    cov_bias_results = np.zeros((len(avals), len(return_vars_covbias)))
    # Loop over A groups
    for k, aval in enumerate(avals):
        A_col = (dat['A'] == aval).astype(int)
        A_col_prob = 1 - dat['Pa'] if aval == 0 else dat['Pa']
        cov_bias_results[k, :] = calculate_bias(dat, A_col=A_col, A_col_prob=A_col_prob, epsilon=epsilon, epsilon_prime=epsilon_prime)
    return cov_bias_results.flatten()

# Calculate bias terms for a specific A=a group
def calculate_bias(dat, A_col, A_col_prob, epsilon, epsilon_prime, p_a_input=None, p_Z_input=None, sens_param_3=None, short_return=False):
    # Placeholder for the actual implementation
    # This function should return a vector including: cov1, cov2, bias (full and denominator), alternate versions of 
    # bias and denominator, true/weighted fnr, true/weighted mu
    pass

# Simulate a population and calculate metrics
def sim_bias_pop(param_1, param_2, param_3, N_pop, p_Y):
    # Placeholder for the actual implementation
    # This function should return a dictionary with 'data' and 'pop_results'
    pass

# Simulate a population and sample to capture both sampling error and bias
def sim_bias_sample(data_pop, N_sample, R, epsilon, epsilon_prime):
    # Placeholder for the actual implementation
    # This function should return a matrix or list depending on the length of epsilon and epsilon_prime
    pass

# Simulate data and get bias terms
def sim_bias(param_1, param_2, param_3, param_4=0, N, p_Y, R, epsilon, epsilon_prime):
    # Placeholder for the actual implementation
    # This function should return a matrix or list depending on the length of epsilon and epsilon_prime
    pass

# Shape simulation results to get means, %-tile intervals
def sim_shape(sim_results, epsilon, epsilon_prime):
    # Aggregate over columns: get means, 95%-tile interval
    mean_vec = np.mean(sim_results, axis=0)
    cilow_vec = np.percentile(sim_results, 2.5, axis=0)
    cihigh_vec = np.percentile(sim_results, 97.5, axis=0)
    res = pd.DataFrame({
        'metric': sim_results.columns,
        'mean': mean_vec,
        'cilow': cilow_vec,
        'cihigh': cihigh_vec
    })
    return res

# Shape simulation results for epsilon subsets
def sim_shape_epsilonsub(sim_results, metrics, param_vary, epsilon_sub, epsilon_prime_sub):
    # Placeholder for the actual implementation
    # This function should return a transformed version of the simulation results for plotting
    pass