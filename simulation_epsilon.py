import numpy as np
import pandas as pd
import pickle
from simulation_functions import sim_bias, sim_shape, sim_shape_epsilonsub
import matplotlib.pyplot as plt

# Load saved BISG accuracy results
with open("rda_files/param_bisg_mean.pkl", "rb") as f:
    param_bisg_mean = pickle.load(f)
with open("rda_files/param_bisg_auc.pkl", "rb") as f:
    param_bisg_auc = pickle.load(f)

#####################
# Setup and run sim
#####################
param1_vals = np.arange(-0.5, 0.6, 0.1)
param2_vals = np.arange(-1, 1.1, 0.1)
param3_vals = np.arange(0, 21, 1)

N = 5000
p_Y = np.array([-0.2])

epsilon_vals = np.arange(-0.05, 0.06, 0.01)
epsilon_prime_vals = np.arange(-0.05, 0.06, 0.01)

# Parameter 1
param1_res = [sim_bias(param_1=p, param_2=0, param_3=20, N=N, p_Y=p_Y, R=100, epsilon=epsilon_vals, epsilon_prime=epsilon_prime_vals) for p in param1_vals]
with open("rda_files/param1_epsilon.pkl", "wb") as f:
    pickle.dump(param1_res, f)

param1_res_long = pd.concat([sim_shape(r, epsilon=epsilon_vals, epsilon_prime=epsilon_prime_vals) for r in param1_res], keys=param1_vals, names=['param1'])
param1_res_long['param1'] = param1_res_long.index.get_level_values('param1').astype(float)
param1_res_long[['metric', 'group']] = param1_res_long['metric'].str.split('.', expand=True)

## Subsets of epsilon_grid
param1_res_sub_02 = sim_shape_epsilonsub(param1_res_long, metrics=["bias_obs_fnr", "bias_epsilon_fnr"], param_vary="param1", epsilon_sub=[-0.02, 0.02], epsilon_prime_sub=[-0.02, 0.02])
param1_res_sub_03 = sim_shape_epsilonsub(param1_res_long, metrics=["bias_obs_fnr", "bias_epsilon_fnr"], param_vary="param1", epsilon_sub=[-0.03, 0.03], epsilon_prime_sub=[-0.03, 0.03])
param1_res_sub_04 = sim_shape_epsilonsub(param1_res_long, metrics=["bias_obs_fnr", "bias_epsilon_fnr"], param_vary="param1", epsilon_sub=[-0.04, 0.04], epsilon_prime_sub=[-0.04, 0.04])

## Transformed version for plot
param1_res_plot = param1_res_long[param1_res_long['metric'].isin(["bias_obs_fnr", "bias_epsilon_fnr"])].groupby(['param1', 'metric', 'group']).agg(
    meanmax=('mean', 'max'),
    meanmin=('mean', 'min'),
    cilowmax=('cilow', 'max'),
    cilowmin=('cilow', 'min'),
    cihighmax=('cihigh', 'max'),
    cihighmin=('cihigh', 'min')
).reset_index().pivot(index=['param1', 'group'], columns='metric', values=['meanmax', 'meanmin', 'cilowmax', 'cilowmin', 'cihighmax', 'cihighmin']).reset_index()

param1_res_plot = param1_res_plot.assign(e_range="05").append([param1_res_sub_02, param1_res_sub_03, param1_res_sub_04], ignore_index=True).assign(
    e_range=lambda df: df['e_range'].map({"05": 0.1, "02": 0.3, "03": 0.5, "04": 0.7})
)

# Plotting
plt.figure()
for group in param1_res_plot['group'].unique():
    subset = param1_res_plot[param1_res_plot['group'] == group]
    plt.fill_between(subset['param1'], subset['meanmin_bias_epsilon_fnr'], subset['meanmax_bias_epsilon_fnr'], alpha=subset['e_range'], label=f'epsilon range {group}')
    plt.plot(subset['param1'], subset['meanmax_bias_obs_fnr'], linestyle='solid', label='Observed')
plt.xlabel('Conditional dependence parameter')
plt.ylabel('FNR')
plt.legend(title='Bias measure')
plt.title('Bias - Conditional Dependence')
plt.show()