import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import multivariate_normal
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

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
        bias_dict = calculate_bias(dat, A_col=A_col, A_col_prob=A_col_prob, epsilon=epsilon, epsilon_prime=epsilon_prime)
        cov_bias_results[k, :] = [bias_dict[var] for var in return_vars_covbias]
    return cov_bias_results.flatten()

def calculate_bias(dat, A_col, A_col_prob, epsilon, epsilon_prime, p_a_input=None, p_Z_input=None, sens_param_3=None, short_return=False):
    # Calculate true and weighted FNR, TPR, FPR, TNR, PPV, NPV
    fnr_true = np.mean((dat['Y'] == 1) & (dat['Yhat'] == 0) * A_col) / np.mean(A_col * dat['Y'])
    fnr_wgt = np.sum((dat['Y'] == 1) & (dat['Yhat'] == 0) * A_col_prob) / np.sum(A_col_prob * dat['Y'])
    
    tpr_true = np.mean((dat['Y'] == 1) & (dat['Yhat'] == 1) * A_col) / np.mean(A_col * dat['Y'])
    tpr_wgt = np.sum((dat['Y'] == 1) & (dat['Yhat'] == 1) * A_col_prob) / np.sum(A_col_prob * dat['Y'])
    
    fpr_true = np.mean((dat['Y'] == 0) & (dat['Yhat'] == 1) * A_col) / np.mean(A_col * (1 - dat['Y']))
    fpr_wgt = np.sum((dat['Y'] == 0) & (dat['Yhat'] == 1) * A_col_prob) / np.sum(A_col_prob * (1 - dat['Y']))
    
    tnr_true = np.mean((dat['Y'] == 0) & (dat['Yhat'] == 0) * A_col) / np.mean(A_col * (1 - dat['Y']))
    tnr_wgt = np.sum((dat['Y'] == 0) & (dat['Yhat'] == 0) * A_col_prob) / np.sum(A_col_prob * (1 - dat['Y']))
    
    ppv_true = np.mean((dat['Yhat'] == 1) & (dat['Y'] == 1) * A_col) / np.mean(A_col * dat['Yhat'])
    ppv_wgt = np.sum((dat['Yhat'] == 1) & (dat['Y'] == 1) * A_col_prob) / np.sum(A_col_prob * dat['Yhat'])
    
    npv_true = np.mean((dat['Yhat'] == 0) & (dat['Y'] == 0) * A_col) / np.mean(A_col * (1 - dat['Yhat']))
    npv_wgt = np.sum((dat['Yhat'] == 0) & (dat['Y'] == 0) * A_col_prob) / np.sum(A_col_prob * (1 - dat['Yhat']))
    
    # Calculate observed bias
    bias_obs_fnr = fnr_wgt - fnr_true
    bias_obs_abs_fnr = np.abs(bias_obs_fnr)
    
    bias_obs_tpr = tpr_wgt - tpr_true
    bias_obs_abs_tpr = np.abs(bias_obs_tpr)
    
    bias_obs_fpr = fpr_wgt - fpr_true
    bias_obs_abs_fpr = np.abs(bias_obs_fpr)
    
    bias_obs_tnr = tnr_wgt - tnr_true
    bias_obs_abs_tnr = np.abs(bias_obs_tnr)
    
    bias_obs_ppv = ppv_wgt - ppv_true
    bias_obs_abs_ppv = np.abs(bias_obs_ppv)
    
    bias_obs_npv = npv_wgt - npv_true
    bias_obs_abs_npv = np.abs(bias_obs_npv)
    
    # Calculate epsilon bias
    fnr_epsilon = fnr_wgt + epsilon
    tpr_epsilon = tpr_wgt + epsilon
    fpr_epsilon = fpr_wgt + epsilon
    tnr_epsilon = tnr_wgt + epsilon
    ppv_epsilon = ppv_wgt + epsilon
    npv_epsilon = npv_wgt + epsilon
    
    # Calculate bias bounds
    bias_bound_abs_fnr = np.abs(fnr_epsilon - fnr_true)
    bias_bound_abs_tpr = np.abs(tpr_epsilon - tpr_true)
    bias_bound_abs_fpr = np.abs(fpr_epsilon - fpr_true)
    bias_bound_abs_tnr = np.abs(tnr_epsilon - tnr_true)
    bias_bound_abs_ppv = np.abs(ppv_epsilon - ppv_true)
    bias_bound_abs_npv = np.abs(npv_epsilon - npv_true)
    
    # Calculate inequality terms
    a1_ineq_leftside_fnr = np.abs(fnr_wgt - fnr_true)
    a1_ineq_rgtside_fnr = np.abs(fnr_epsilon - fnr_true)
    
    a1_ineq_leftside_tpr = np.abs(tpr_wgt - tpr_true)
    a1_ineq_rgtside_tpr = np.abs(tpr_epsilon - tpr_true)
    
    a1_ineq_leftside_fpr = np.abs(fpr_wgt - fpr_true)
    a1_ineq_rgtside_fpr = np.abs(fpr_epsilon - fpr_true)
    
    a1_ineq_leftside_tnr = np.abs(tnr_wgt - tnr_true)
    a1_ineq_rgtside_tnr = np.abs(tnr_epsilon - tnr_true)
    
    a1_ineq_leftside_ppv = np.abs(ppv_wgt - ppv_true)
    a1_ineq_rgtside_ppv = np.abs(ppv_epsilon - ppv_true)
    
    a1_ineq_leftside_npv = np.abs(npv_wgt - npv_true)
    a1_ineq_rgtside_npv = np.abs(npv_epsilon - npv_true)
    
    # Calculate true bias
    bias_true_fnr = fnr_wgt - fnr_true
    
    # Combine results into a dictionary
    results = {
        'fnr_true': fnr_true, 'fnr_wgt': fnr_wgt, 'tpr_true': tpr_true, 'tpr_wgt': tpr_wgt,
        'fpr_true': fpr_true, 'fpr_wgt': fpr_wgt, 'tnr_true': tnr_true, 'tnr_wgt': tnr_wgt,
        'ppv_true': ppv_true, 'ppv_wgt': ppv_wgt, 'npv_true': npv_true, 'npv_wgt': npv_wgt,
        'bias_obs_fnr': bias_obs_fnr, 'bias_obs_abs_fnr': bias_obs_abs_fnr, 'bias_obs_tpr': bias_obs_tpr,
        'bias_obs_abs_tpr': bias_obs_abs_tpr, 'bias_obs_fpr': bias_obs_fpr, 'bias_obs_abs_fpr': bias_obs_abs_fpr,
        'bias_obs_tnr': bias_obs_tnr, 'bias_obs_abs_tnr': bias_obs_abs_tnr, 'bias_obs_ppv': bias_obs_ppv,
        'bias_obs_abs_ppv': bias_obs_abs_ppv, 'bias_obs_npv': bias_obs_npv, 'bias_obs_abs_npv': bias_obs_abs_npv,
        'fnr_epsilon': fnr_epsilon, 'tpr_epsilon': tpr_epsilon, 'fpr_epsilon': fpr_epsilon, 'tnr_epsilon': tnr_epsilon,
        'ppv_epsilon': ppv_epsilon, 'npv_epsilon': npv_epsilon, 'bias_true_fnr': bias_true_fnr,
        'bias_bound_abs_fnr': bias_bound_abs_fnr, 'bias_bound_abs_tpr': bias_bound_abs_tpr, 'bias_bound_abs_fpr': bias_bound_abs_fpr,
        'bias_bound_abs_tnr': bias_bound_abs_tnr, 'bias_bound_abs_ppv': bias_bound_abs_ppv, 'bias_bound_abs_npv': bias_bound_abs_npv,
        'a1_ineq_leftside_fnr': a1_ineq_leftside_fnr, 'a1_ineq_rgtside_fnr': a1_ineq_rgtside_fnr,
        'a1_ineq_leftside_tpr': a1_ineq_leftside_tpr, 'a1_ineq_rgtside_tpr': a1_ineq_rgtside_tpr,
        'a1_ineq_leftside_fpr': a1_ineq_leftside_fpr, 'a1_ineq_rgtside_fpr': a1_ineq_rgtside_fpr,
        'a1_ineq_leftside_tnr': a1_ineq_leftside_tnr, 'a1_ineq_rgtside_tnr': a1_ineq_rgtside_tnr,
        'a1_ineq_leftside_ppv': a1_ineq_leftside_ppv, 'a1_ineq_rgtside_ppv': a1_ineq_rgtside_ppv,
        'a1_ineq_leftside_npv': a1_ineq_leftside_npv, 'a1_ineq_rgtside_npv': a1_ineq_rgtside_npv,
        'bias_epsilon_fnr': bias_bound_abs_fnr, 'bias_epsilon_tpr': bias_bound_abs_tpr, 'bias_epsilon_fpr': bias_bound_abs_fpr,
        'bias_epsilon_tnr': bias_bound_abs_tnr, 'bias_epsilon_ppv': bias_bound_abs_ppv, 'bias_epsilon_npv': bias_bound_abs_npv
    }
    
    return results
# Simulate a population and calculate metrics
def sim_bias_pop(param_1, param_2, param_3, N_pop, p_Y):
    # Simulate Z
    Z = np.random.normal(-0.4, 1, N_pop)
    
    # Simulate A and P(A=a|Z)
    Sigma_mat = np.array([[20, param_3], [param_3, 20]])
    pA_mat = np.array([multivariate_normal.rvs(mean=[z, z + param_2], cov=Sigma_mat) for z in Z])
    A = np.random.binomial(1, expit(pA_mat[:, 0]))
    
    dat = pd.DataFrame({
        'A': A,
        'Pa': expit(pA_mat[:, 1]),
        'Z': Z
    })
    
    # Simulate X: one component dependent on Z
    dat['X_2'] = np.random.normal(0, 0.5, N_pop)
    dat['X_3'] = np.random.normal(1, 0.5, N_pop)
    dat['X_4'] = np.random.normal(-1, 0.5, N_pop)
    
    # Simulate Y
    prob_Y = expit(p_Y + Z + dat['X_2'] + dat['X_3'] + dat['X_4'] + param_1 * A)
    dat['Y'] = np.random.binomial(1, prob_Y)
    
    # Split training and test data
    dat_train = dat.sample(frac=0.25, random_state=1)
    dat_test = dat.drop(dat_train.index)
    
    # Yhat based on model trained on training data
    model_Y = LogisticRegression()
    model_Y.fit(dat_train[['A', 'Z', 'X_2', 'X_3', 'X_4']], dat_train['Y'])
    Yhat_prob = model_Y.predict_proba(dat_test[['A', 'Z', 'X_2', 'X_3', 'X_4']])[:, 1]
    dat_test['Yhat'] = (Yhat_prob >= 0.5).astype(int)
    dat = dat_test
    
    # Get metrics on whole population
    pop_results = get_avals_bias(dat, epsilon=0.05, epsilon_prime=0.05)
    pop_results = {f"{var}.{suffix}": val for var, val in zip(return_vars_covbias, pop_results) for suffix in ["0", "1"]}
    
    return {'data': dat, 'pop_results': pop_results}

# Simulate a population and sample to capture both sampling error and bias
def sim_bias_sample(data_pop, N_sample, R, epsilon, epsilon_prime):
    N_pop = len(data_pop)
    
    # Results matrix/list for samples
    if len(epsilon) == 1 and len(epsilon_prime) == 1:
        samp_results = np.zeros((R, len(return_vars_covbias) * 2))
        colnames = [f"{var}.{suffix}" for var in return_vars_covbias for suffix in ["0", "1"]]
    else:
        epsilon_grid = pd.DataFrame([(e, ep) for e in epsilon for ep in epsilon_prime], columns=['epsilon', 'epsilon_prime'])
        samp_results = [None] * R
    
    # Draw samples and calculate metrics
    for i in range(R):
        sample_ind = np.random.choice(N_pop, N_sample, replace=False)
        dat_samp = data_pop.iloc[sample_ind]
        
        # Get bias terms
        if len(epsilon) == 1 and len(epsilon_prime) == 1:
            bias_res = get_avals_bias(dat_samp, epsilon=epsilon[0], epsilon_prime=epsilon_prime[0])
            samp_results[i, :] = bias_res
        else:
            e_mat = np.zeros((len(epsilon_grid), len(return_vars_covbias) * 2))
            colnames = [f"{var}.{suffix}" for var in return_vars_covbias for suffix in ["0", "1"]]
            # Loop over epsilon grid values
            for j, (e, e_prime) in enumerate(epsilon_grid.values):
                bias_res = get_avals_bias(dat_samp, epsilon=e, epsilon_prime=e_prime)
                e_mat[j, :] = bias_res
            # Add matrix of epsilon, epsilon' results to replications list
            samp_results[i] = e_mat
    
    return samp_results

# Simulate data and get bias terms
def sim_bias(param_1, param_2, param_3, param_4=0, N=1000, p_Y=-0.2, R=100, epsilon=[0.05], epsilon_prime=[0.05]):
    # Results matrix/list
    if len(epsilon) == 1 and len(epsilon_prime) == 1:
        sim_results = np.zeros((R, len(return_vars_covbias) * 2))
        colnames = [f"{var}.{suffix}" for var in return_vars_covbias for suffix in ["0", "1"]]
    else:
        epsilon_grid = pd.DataFrame([(e, ep) for e in epsilon for ep in epsilon_prime], columns=['epsilon', 'epsilon_prime'])
        sim_results = [None] * R
    
    # Vector for AUC values
    auc_results = np.zeros(R)
    
    for i in range(R):
        # Simulate Z
        Z = np.random.normal(-0.4, 1, N)
        
        # Simulate A and P(A=a|Z)
        Sigma_mat = np.array([[20, param_3], [param_3, 20]])
        pA_mat = np.array([multivariate_normal.rvs(mean=[z, z + param_2], cov=Sigma_mat) for z in Z])
        A = np.random.binomial(1, expit(pA_mat[:, 0]))
        
        dat = pd.DataFrame({
            'A': A,
            'Pa': expit(pA_mat[:, 1]),
            'Z': Z
        })
        
        # Calculate AUC
        auc_results[i] = roc_auc_score(dat['A'], dat['Pa'])
        
        # Simulate X: one component dependent on Z
        dat['X_2'] = np.random.normal(0, 0.5, N)
        dat['X_3'] = np.random.normal(1, 0.5, N)
        dat['X_4'] = np.random.normal(-1, 0.5, N)
        
        # Simulate Y
        prob_Y = expit(p_Y + Z + dat['X_2'] + dat['X_3'] + dat['X_4'] + param_1 * A)
        dat['Y'] = np.random.binomial(1, prob_Y)
        
        # Split training and test data
        dat_train = dat.sample(frac=0.25, random_state=1)
        dat_test = dat.drop(dat_train.index)
        
        # Yhat based on model trained on training data
        model_Y = LogisticRegression()
        model_Y.fit(dat_train[['A', 'Z', 'X_2', 'X_3', 'X_4']], dat_train['Y'])
        Yhat_prob = model_Y.predict_proba(dat_test[['A', 'Z', 'X_2', 'X_3', 'X_4']])[:, 1]
        Yhat_prob += np.random.normal(0, param_4, len(Yhat_prob))
        dat_test['Yhat'] = (Yhat_prob >= 0.5).astype(int)
        dat = dat_test
        
        # Get bias terms
        if len(epsilon) == 1 and len(epsilon_prime) == 1:
            bias_res = get_avals_bias(dat, epsilon=epsilon[0], epsilon_prime=epsilon_prime[0])
            sim_results[i, :] = bias_res
        else:
            e_mat = np.zeros((len(epsilon_grid), len(return_vars_covbias) * 2))
            colnames = [f"{var}.{suffix}" for var in return_vars_covbias for suffix in ["0", "1"]]
            # Loop over epsilon grid values
            for j, (e, e_prime) in enumerate(epsilon_grid.values):
                bias_res = get_avals_bias(dat, epsilon=e, epsilon_prime=e_prime)
                e_mat[j, :] = bias_res
            # Add matrix of epsilon, epsilon' results to replications list
            sim_results[i] = e_mat
    
    return {'sim_results': sim_results, 'auc_results': auc_results}

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
def sim_shape_epsilonsub(sim_results_long, metrics, param_vary, epsilon_sub, epsilon_prime_sub):
    # Round epsilon and epsilon_prime columns
    sim_results_long['e_col'] = sim_results_long['e_col'].round(2)
    sim_results_long['eprime_col'] = sim_results_long['eprime_col'].round(2)
    
    # Filter based on epsilon and epsilon_prime ranges
    filtered_results = sim_results_long[
        (sim_results_long['e_col'] >= min(epsilon_sub)) & (sim_results_long['e_col'] <= max(epsilon_sub)) &
        (sim_results_long['eprime_col'] >= min(epsilon_prime_sub)) & (sim_results_long['eprime_col'] <= max(epsilon_prime_sub)) &
        (sim_results_long['metric'].isin(metrics))
    ]
    
    # Group by param_vary, metric, and group, then summarize
    grouped_results = filtered_results.groupby([param_vary, 'metric', 'group']).agg(
        meanmax=('mean', 'max'),
        meanmin=('mean', 'min'),
        cilowmax=('cilow', 'max'),
        cilowmin=('cilow', 'min'),
        cihighmax=('cihigh', 'max'),
        cihighmin=('cihigh', 'min')
    ).reset_index()
    
    # Pivot wider
    pivoted_results = grouped_results.pivot_table(
        index=[param_vary, 'group'],
        columns='metric',
        values=['meanmax', 'meanmin', 'cilowmax', 'cilowmin', 'cihighmax', 'cihighmin'],
        aggfunc='first'
    ).reset_index()
    
    # Flatten the multi-level columns
    pivoted_results.columns = ['_'.join(col).strip() if col[1] else col[0] for col in pivoted_results.columns.values]
    
    # Add e_range column
    pivoted_results['e_range'] = str(max(epsilon_sub)).replace("0.", "")
    
    return pivoted_results