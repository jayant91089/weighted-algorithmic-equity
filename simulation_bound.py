import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.stats import norm
import pickle

# Load saved BISG accuracy results
with open("rda_files/param_bisg_mean.pkl", "rb") as f:
    param_bisg_mean = pickle.load(f)
with open("rda_files/param_bisg_auc.pkl", "rb") as f:
    param_bisg_auc = pickle.load(f)

#####################
# Setup and run sim
#####################
facet_labs_A = {0: "A=0", 1: "A=1"}

param1_vals = np.arange(-0.5, 0.6, 0.1)
param2_vals = np.arange(-1, 1.1, 0.1)
param3_vals = np.arange(0, 21, 1)

N = 5000
p_Y = np.array([-0.2])

epsilon_vals = np.arange(-0.05, 0.06, 0.01)
epsilon_prime_vals = np.arange(-0.05, 0.06, 0.01)

# Load saved results
with open("rda_files/param1_epsilon.pkl", "rb") as f:
    param1_res = pickle.load(f)
with open("rda_files/param2_epsilon.pkl", "rb") as f:
    param2_res = pickle.load(f)
with open("rda_files/param3_epsilon.pkl", "rb") as f:
    param3_res = pickle.load(f)

# Shape results
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

param1_res_long = pd.concat([sim_shape(r, epsilon_vals, epsilon_prime_vals) for r in param1_res], keys=param1_vals, names=['param1'])
param1_res_long['param1'] = param1_res_long.index.get_level_values('param1').astype(float)
param1_res_long[['metric', 'group']] = param1_res_long['metric'].str.split('.', expand=True)

param2_res_long = pd.concat([sim_shape(r, epsilon_vals, epsilon_prime_vals) for r in param2_res], keys=param2_vals, names=['param2'])
param2_res_long['param2'] = param2_res_long.index.get_level_values('param2').astype(float)
param2_res_long[['metric', 'group']] = param2_res_long['metric'].str.split('.', expand=True)

param3_res_long = pd.concat([sim_shape(r, epsilon_vals, epsilon_prime_vals) for r in param3_res], keys=param3_vals, names=['param3'])
param3_res_long['param3'] = param3_res_long.index.get_level_values('param3').astype(float)
param3_res_long[['metric', 'group']] = param3_res_long['metric'].str.split('.', expand=True)

test_auc_3 = [r['auc_results'][0] for r in param3_res]

# Transformed versions for plots
def transform_results(res_long, metrics):
    res_plot = res_long[res_long['metric'].isin(metrics)].groupby(['param1', 'metric', 'group']).agg(
        meanmax=('mean', 'max'),
        meanmin=('mean', 'min'),
        cilowmax=('cilow', 'max'),
        cilowmin=('cilow', 'min'),
        cihighmax=('cihigh', 'max'),
        cihighmin=('cihigh', 'min')
    ).reset_index().pivot(index=['param1', 'group'], columns='metric', values=['meanmax', 'meanmin', 'cilowmax', 'cilowmin', 'cihighmax', 'cihighmin']).reset_index()
    return res_plot

param1_res_plot = transform_results(param1_res_long, ["bias_bound_abs_fnr", "bias_obs_abs_fnr", "bias_bound_abs_fpr", "bias_obs_abs_fpr",
                                                      "bias_bound_abs_ppv", "bias_obs_abs_ppv", "bias_bound_abs_npv", "bias_obs_abs_npv",
                                                      "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr",
                                                      "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr",
                                                      "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                                                      "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv"])

param2_res_plot = transform_results(param2_res_long, ["bias_bound_abs_fnr", "bias_obs_abs_fnr", "bias_bound_abs_fpr", "bias_obs_abs_fpr",
                                                      "bias_bound_abs_ppv", "bias_obs_abs_ppv", "bias_bound_abs_npv", "bias_obs_abs_npv",
                                                      "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr", "pa_ratio",
                                                      "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr",
                                                      "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                                                      "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv"])

param3_res_plot = transform_results(param3_res_long, ["bias_bound_abs_fnr", "bias_obs_abs_fnr", "bias_bound_abs_fpr", "bias_obs_abs_fpr",
                                                      "bias_bound_abs_ppv", "bias_obs_abs_ppv", "bias_bound_abs_npv", "bias_obs_abs_npv",
                                                      "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr",
                                                      "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr",
                                                      "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                                                      "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv"])

# Combine parameters 2 and 3
param23_list = {'param2': param2_res_plot.rename(columns={'param2': 'param'}), 'param3': param3_res_plot.rename(columns={'param3': 'param'})}
param23_res = pd.concat(param23_list, names=['param_set'])

# Look up BISG AUC values and get corresponding parameter 3
param3_auc_df = pd.DataFrame({'param3': param3_vals, 'auc': test_auc_3})
param_cov = param3_auc_df['param3'][np.searchsorted(param3_auc_df['auc'], param_bisg_auc)]

param_bisg = pd.DataFrame({
    'param_set': ['param2', 'param3'],
    'param_min': [param_bisg_mean.min(), param_cov.min()],
    'param_max': [param_bisg_mean.max(), param_cov.max()]
})

########
# Plots
########
# Part 1: Conditional dependence
## Bias bound
plt.figure()
for group in param1_res_plot['group'].unique():
    subset = param1_res_plot[param1_res_plot['group'] == group]
    plt.plot(subset['param1'], subset['meanmax_bias_bound_abs_fnr'], linestyle='dashed', label='Proposed bound')
    plt.plot(subset['param1'], subset['meanmax_bias_obs_abs_fnr'], linestyle='solid', label='Observed')
plt.xlabel('Conditional dependence parameter')
plt.ylabel('Bias')
plt.legend(title='Bias measure')
plt.title('Bias Bound - Conditional Dependence')
plt.show()

## Assumption
plt.figure()
for group in param1_res_plot['group'].unique():
    subset = param1_res_plot[param1_res_plot['group'] == group]
    plt.plot(subset['param1'], subset['meanmax_a1_ineq_leftside_fnr'], linestyle='solid', label='|δ|')
    plt.plot(subset['param1'], subset['meanmax_a1_ineq_rgtside_fnr'], linestyle='dashed', label='|δ*|')
plt.xlabel('Conditional dependence parameter')
plt.ylabel('Mean of 100 replications')
plt.legend(title='')
plt.title('Assumption - Conditional Dependence')
plt.show()

## Other metrics
def plot_metric(param_res_plot, metric, ylabel):
    plt.figure()
    for group in param_res_plot['group'].unique():
        subset = param_res_plot[param_res_plot['group'] == group]
        plt.plot(subset['param1'], subset[f'meanmax_{metric}'], linestyle='dashed', label='Proposed bound')
        plt.plot(subset['param1'], subset[f'meanmax_bias_obs_abs_{metric}'], linestyle='solid', label='Observed')
    plt.xlabel('Conditional dependence parameter')
    plt.ylabel(ylabel)
    plt.legend(title='Bias measure')
    plt.title(f'{ylabel} - Conditional Dependence')
    plt.show()

plot_metric(param1_res_plot, 'bias_bound_abs_fpr', 'FPR Bias')
plot_metric(param1_res_plot, 'bias_bound_abs_ppv', 'PPV Bias')
plot_metric(param1_res_plot, 'bias_bound_abs_npv', 'NPV Bias')

# Part 2: Probability ratio and AUC
## Probability ratio
plt.figure()
for group in param2_res_plot['group'].unique():
    subset = param2_res_plot[param2_res_plot['group'] == group]
    plt.plot(subset['param2'], subset['meanmax_pa_ratio'], label=group)
plt.xlabel('p_A accuracy parameter: mean')
plt.ylabel('Ratio')
plt.legend(title='Group')
plt.title('Probability Ratio')
plt.show()

## AUC
plt.figure()
plt.scatter(param3_vals, test_auc_3)
plt.xlabel('p_A accuracy parameter: covariance')
plt.ylabel('AUC')
plt.title('AUC')
plt.show()

## Bias bound
plt.figure()
for param_set in param23_res['param_set'].unique():
    subset = param23_res[param23_res['param_set'] == param_set]
    plt.plot(subset['param'], subset['meanmax_bias_bound_abs_fnr'], linestyle='dashed', label='Proposed bound')
    plt.plot(subset['param'], subset['meanmax_bias_obs_abs_fnr'], linestyle='solid', label='Observed')
    plt.axvline(x=param_bisg[param_bisg['param_set'] == param_set]['param_min'].values[0], color='red')
    plt.axvline(x=param_bisg[param_bisg['param_set'] == param_set]['param_max'].values[0], color='red')
plt.xlabel('p_A accuracy parameter')
plt.ylabel('Bias')
plt.legend(title='Bias measure')
plt.title('Bias Bound - Probability Ratio and AUC')
plt.show()

## Assumption
plt.figure()
for param_set in param23_res['param_set'].unique():
    subset = param23_res[param23_res['param_set'] == param_set]
    plt.plot(subset['param'], subset['meanmax_a1_ineq_leftside_fnr'], linestyle='solid', label='|δ|')
    plt.plot(subset['param'], subset['meanmax_a1_ineq_rgtside_fnr'], linestyle='dashed', label='|δ*|')
plt.xlabel('p_A accuracy parameter: mean')
plt.ylabel('Mean of 100 replications')
plt.legend(title='')
plt.title('Assumption - Probability Ratio and AUC')
plt.show()

## Other metrics
plot_metric(param23_res, 'bias_bound_abs_fpr', 'FPR Bias')
plot_metric(param23_res, 'bias_bound_abs_ppv', 'PPV Bias')
plot_metric(param23_res, 'bias_bound_abs_npv', 'NPV Bias')