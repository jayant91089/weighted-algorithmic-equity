import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

########
# Setup
########
fnr_wgt = np.array([0.35, 0.48])
mean_fn_overall = 0.16
mean_pos_overall = 0.42
fnr_overall = mean_fn_overall / mean_pos_overall
mean_tp_overall = mean_pos_overall - mean_fn_overall
mean_est_fn = np.array([0.09, 0.07])
mean_est_pos = mean_est_fn / fnr_wgt
mean_est_tp = mean_est_pos - mean_est_fn
prop_A = np.array([0.48, 0.52])

# Unknown parameter grid
mean_true_fn_vec = np.arange(0.005, 1.005, 0.005)
mean_true_pos_vec = np.arange(0.005, 1.005, 0.005)
mean_true_grid = pd.DataFrame([(fn, pos) for fn in mean_true_fn_vec for pos in mean_true_pos_vec if fn <= pos], columns=['mean_true_fn_vec', 'mean_true_pos_vec'])

bias_mat = np.zeros((len(mean_true_grid), 2))

group_colors = ["#440154", "#2A788E", "#22A884", "#7AD151"]
group_labels = {1: "Group 1", 2: "Group 2"}

#################################################
# Loop over groups and unknown parameter values
#################################################
for i in range(2):
    for j in range(len(mean_true_grid)):
        ## Set unknown parameters
        mean_true_fn = mean_true_grid.iloc[j, 0]
        mean_true_pos = mean_true_grid.iloc[j, 1]
        ## Calculate bias
        bias_val = (mean_est_fn[i] - mean_true_fn) / mean_true_pos + fnr_wgt[i] * (1 - mean_est_pos[i] / mean_true_pos)
        
        bias_mat[j, i] = bias_val

bias_mat_red = pd.DataFrame(bias_mat, columns=["G_1", "G_2"])

# Remove impossible values (higher than overall group proportions OR mean_true_fn > mean_true_pos)
for i in range(2):
    temp_df = mean_true_grid.copy()
    temp_df['is_impossible'] = (temp_df['mean_true_fn_vec'] > mean_fn_overall) | (temp_df['mean_true_pos_vec'] > mean_pos_overall) | (temp_df['mean_true_fn_vec'] > temp_df['mean_true_pos_vec'])
    temp_df['epsilon'] = (mean_est_fn[i] - temp_df['mean_true_fn_vec']) / mean_fn_overall
    temp_df['epsilon_prime'] = (mean_est_tp[i] - (temp_df['mean_true_pos_vec'] - temp_df['mean_true_fn_vec'])) / mean_tp_overall

    # Filter out impossible values
    temp_df = temp_df[~temp_df['is_impossible']]

    # Plotting
    plt.figure()
    plt.scatter(temp_df['epsilon'], temp_df['epsilon_prime'], c=bias_mat_red.iloc[temp_df.index, i], cmap='viridis')
    plt.colorbar(label='Bias')
    plt.xlabel('Epsilon')
    plt.ylabel('Epsilon Prime')
    plt.title(f'Bias for Group {i+1}')
    plt.show()