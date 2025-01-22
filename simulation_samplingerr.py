import numpy as np
import pandas as pd
import pickle
from simulation_functions import sim_bias_pop, sim_bias_sample, sim_shape
import matplotlib.pyplot as plt

#####################
# Setup and run sim
#####################
N_vals = np.arange(1000, 10001, 1000)
param1_vals = np.arange(-0.5, 0.6, 0.1)

p_Y = np.array([-0.2])

epsilon_vals = np.arange(-0.05, 0.06, 0.01)
epsilon_prime_vals = np.arange(-0.05, 0.06, 0.01)

################ Sample size
simN_res_pop = sim_bias_pop(param_1=0.25, param_2=0, param_3=20, N_pop=50000, p_Y=p_Y)

simN_res_all = [sim_bias_sample(data_pop=simN_res_pop['data'], N_sample=n, R=100, epsilon=epsilon_vals, epsilon_prime=epsilon_prime_vals) for n in N_vals]

simN_res_all = [{'samp_results': r, 'pop_results': simN_res_pop['pop_results']} for r in simN_res_all]

simN_res_long = pd.concat([sim_shape(r, epsilon=epsilon_vals, epsilon_prime=epsilon_prime_vals) for r in simN_res_all], keys=N_vals, names=['N'])
simN_res_long['N'] = simN_res_long.index.get_level_values('N').astype(float)
simN_res_long[['metric', 'group']] = simN_res_long['metric'].str.split('.', expand=True)

## Transformed version for plot
simN_bias_plot = simN_res_long[(simN_res_long['metric'] == 'bias_obs_fnr') & (simN_res_long['e_col'] == 0.05) & (simN_res_long['eprime_col'] == 0.05)].drop(columns=['pop', 'e_col', 'eprime_col'])

################ Conditional dependence: N_sample = 2,000
simparam1_res_all_2K = []
for p in param1_vals:
    res_pop = sim_bias_pop(param_1=p, param_2=0, param_3=20, N_pop=50000, p_Y=p_Y)
    res_samp = sim_bias_sample(data_pop=res_pop['data'], N_sample=2000, R=100, epsilon=epsilon_vals, epsilon_prime=epsilon_prime_vals)
    simparam1_res_all_2K.append({'samp_results': res_samp, 'pop_results': res_pop['pop_results']})

simparam1_res_long_2K = pd.concat([sim_shape(r, epsilon=epsilon_vals, epsilon_prime=epsilon_prime_vals) for r in simparam1_res_all_2K], keys=param1_vals, names=['param1'])
simparam1_res_long_2K['param1'] = simparam1_res_long_2K.index.get_level_values('param1').astype(float)
simparam1_res_long_2K[['metric', 'group']] = simparam1_res_long_2K['metric'].str.split('.', expand=True)

## Transformed version for plot
simparam1_bias_plot_2K = simparam1_res_long_2K[(simparam1_res_long_2K['metric'] == 'bias_obs_fnr') & (simparam1_res_long_2K['e_col'] == 0.05) & (simparam1_res_long_2K['eprime_col'] == 0.05)].drop(columns=['pop', 'e_col', 'eprime_col'])

#########
# Plots
#########
# Sample size
plt.figure()
for group in simN_bias_plot['group'].unique():
    subset = simN_bias_plot[simN_bias_plot['group'] == group]
    plt.errorbar(subset['N'], subset['mean'], yerr=[subset['mean'] - subset['cilow'], subset['cihigh'] - subset['mean']], fmt='o', label=group)
plt.axhline(y=0, color='black', linestyle='solid')
plt.xlabel('Sample size')
plt.ylabel('Bias')
plt.legend(title='Group')
plt.title('Bias by Sample Size')
plt.show()

# Conditional dependence
plt.figure()
for group in simparam1_bias_plot_2K['group'].unique():
    subset = simparam1_bias_plot_2K[simparam1_bias_plot_2K['group'] == group]
    plt.errorbar(subset['param1'], subset['mean'], yerr=[subset['mean'] - subset['cilow'], subset['cihigh'] - subset['mean']], fmt='o', label=group)
plt.axhline(y=0, color='black', linestyle='solid')
plt.xlabel('Conditional dependence parameter')
plt.ylabel('Bias')
plt.legend(title='Group')
plt.title('Bias by Conditional Dependence Parameter')
plt.show()