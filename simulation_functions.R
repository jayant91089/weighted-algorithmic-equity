library(tidyverse)
expit <- function(x){exp(x)/(1+exp(x))}
p_trunc <- function(probs) {
  pmin(pmax(probs, 0.005), 0.995)
}

# Variables returned by calc_cov_bias
return_vars_covbias <- c("fnr_true", "fnr_wgt", "tpr_true", "tpr_wgt", "fpr_true", "fpr_wgt", "tnr_true", "tnr_wgt", "ppv_true", "ppv_wgt", "npv_true", "npv_wgt",
                         "bias_obs_fnr", "bias_obs_abs_fnr", "bias_obs_tpr", "bias_obs_abs_tpr", "bias_obs_fpr", "bias_obs_abs_fpr", "bias_obs_tnr", "bias_obs_abs_tnr", 
                         "bias_obs_ppv", "bias_obs_abs_ppv", "bias_obs_npv", "bias_obs_abs_npv",
                         "fnr_epsilon", "tpr_epsilon", "fpr_epsilon", "tnr_epsilon", "ppv_epsilon", "npv_epsilon",
                         "pa_ratio", "bias_true_fnr",
                         "bias_bound_abs_fnr", "bias_bound_abs_tpr", "bias_bound_abs_fpr", "bias_bound_abs_tnr", "bias_bound_abs_ppv", "bias_bound_abs_npv",
                         "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr", "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                         "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr", "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv",
                         "bias_epsilon_fnr", "bias_epsilon_tpr", "bias_epsilon_fpr", "bias_epsilon_tnr", "bias_epsilon_ppv", "bias_epsilon_npv")

# Helper function to get bias terms across multiple A groups ##############
## Inputs:
### data: input data set
### epsilon, epsilon': sensitivity parameters for estimating bias

## Outputs: vector of bias terms including all A groups (in factor level order)
#############################################################################
get_avals_bias <- function(dat, epsilon, epsilon_prime) {
  avals <- as.numeric(levels(as.factor(dat$A)))
  cov_bias_results <- matrix(nrow = length(avals), ncol = length(return_vars_covbias))
  # Loop over A groups
  for(k in 1:length(avals)) {
    A_col <- as.numeric(dat$A == avals[k])
    if(avals[k]==0) { A_col_prob <- 1-dat$Pa } else { A_col_prob <- dat$Pa }
    cov_bias_results[k,] <- calculate_bias(dat, A_col=A_col, A_col_prob=A_col_prob, epsilon=epsilon, epsilon_prime=epsilon_prime)
  }
  return(as.vector(t(cov_bias_results)))
}

# Calculate bias terms for a specific A=a group #################################
## Inputs:
### dat: simulated data with 
### A_col: indicator variable for A=a
### A_col_prob: vector of probabilities for P(A=a)
### p_a_input: P(A=a|Z=z)
### p_Z_input: P(Z=z)
### epsilon, epsilon': sensitivity parameters for estimating bias
### sens_param_3: third sensitivity parameter
### short_return: (T/F) use shortened return vector to accommodate looping over many A categories

## Outputs:
### named vector including: cov1, cov2, bias (full and denominator), alternate versions of 
#### bias and denominator, true/weighted fnr, true/weighted mu
####################################################################################################################
calculate_bias <- function(dat, A_col, A_col_prob, epsilon, epsilon_prime) {
  # Performance metrics
  ## FNR
  fnr_true <- mean(dat$Y*(1-dat$Yhat)*A_col)/mean(A_col*dat$Y)
  fnr_wgt <- mean(dat$Y*(1-dat$Yhat)*A_col_prob)/mean(A_col_prob*dat$Y)
  fnr_marg <- mean(dat$Y*(1-dat$Yhat))/mean(dat$Y)
  ## TPR (1-FNR)
  tpr_true <- mean(dat$Y*dat$Yhat*A_col)/mean(A_col*dat$Y)
  tpr_wgt <- mean(dat$Y*dat$Yhat*A_col_prob)/mean(A_col_prob*dat$Y)
  tpr_marg <- mean(dat$Y*dat$Yhat)/mean(dat$Y)
  ## FPR
  fpr_true <- mean((1-dat$Y)*dat$Yhat*A_col)/mean(A_col*(1-dat$Y))
  fpr_wgt <- mean((1-dat$Y)*dat$Yhat*A_col_prob)/mean(A_col_prob*(1-dat$Y))
  fpr_marg <- mean((1-dat$Y)*dat$Yhat)/mean((1-dat$Y))
  ## TNR (1-FPR)
  tnr_true <- mean((1-dat$Y)*(1-dat$Yhat)*A_col)/mean(A_col*(1-dat$Y))
  tnr_wgt <- mean((1-dat$Y)*(1-dat$Yhat)*A_col_prob)/mean(A_col_prob*(1-dat$Y))
  tnr_marg <- mean((1-dat$Y)*(1-dat$Yhat))/mean((1-dat$Y))
  ## PPV
  ppv_true <- mean(dat$Y*dat$Yhat*A_col)/mean(A_col*dat$Yhat)
  ppv_wgt <- mean(dat$Y*dat$Yhat*A_col_prob)/mean(A_col_prob*dat$Yhat)
  ppv_marg <- mean(dat$Y*dat$Yhat)/mean(dat$Yhat)
  ## NPV
  npv_true <- mean((1-dat$Y)*(1-dat$Yhat)*A_col)/mean(A_col*(1-dat$Yhat))
  npv_wgt <- mean((1-dat$Y)*(1-dat$Yhat)*A_col_prob)/mean(A_col_prob*(1-dat$Yhat))
  npv_marg <- mean((1-dat$Y)*(1-dat$Yhat))/mean(1-dat$Yhat)
   
  # Observed bias
  bias_obs_fnr <- fnr_wgt - fnr_true
  bias_obs_abs_fnr <- abs(bias_obs_fnr)
  bias_obs_fpr <- fpr_wgt - fpr_true
  bias_obs_abs_fpr <- abs(bias_obs_fpr)
  bias_obs_tpr <- tpr_wgt - tpr_true
  bias_obs_abs_tpr <- abs(bias_obs_tpr)
  bias_obs_tnr <- tnr_wgt - tnr_true
  bias_obs_abs_tnr <- abs(bias_obs_tnr)
  
  bias_obs_ppv <- ppv_wgt - ppv_true
  bias_obs_abs_ppv <- abs(bias_obs_ppv)
  bias_obs_npv <- npv_wgt - npv_true
  bias_obs_abs_npv <- abs(bias_obs_npv)
  
  # Probability ratio
  pa_ratio <- mean(dat$Y*A_col_prob)/mean(dat$Y*A_col)
   
  # Analytic expression for bias
  bias_true_fnr <- (mean(A_col_prob*(1-dat$Yhat)*dat$Y) - mean(A_col*(1-dat$Yhat)*dat$Y))/mean(A_col*dat$Y) + fnr_wgt*(1-pa_ratio)
  
  # Bias bound
  bias_bound_abs_fnr <- (1+fnr_wgt)*abs(1- mean(dat$Y*A_col_prob)/mean(dat$Y*A_col))
  bias_bound_abs_tpr <- (1+tpr_wgt)*abs(1- mean(dat$Y*A_col_prob)/mean(dat$Y*A_col))
  
  bias_bound_abs_fpr <- (1+fpr_wgt)*abs(1- mean((1-dat$Y)*A_col_prob)/mean((1-dat$Y)*A_col))
  bias_bound_abs_tnr <- (1+tnr_wgt)*abs(1- mean((1-dat$Y)*A_col_prob)/mean((1-dat$Y)*A_col))
  
  bias_bound_abs_ppv <- (1+ppv_wgt)*abs(1- mean(dat$Yhat*A_col_prob)/mean(dat$Yhat*A_col))
  bias_bound_abs_npv <- (1+npv_wgt)*abs(1- mean((1-dat$Yhat)*A_col_prob)/mean((1-dat$Yhat)*A_col))
 
  ## Terms to check assumption
  a1_ineq_leftside_fnr <- abs(mean(dat$Y*(1-dat$Yhat)*A_col_prob) - mean(dat$Y*(1-dat$Yhat)*A_col))
  a1_ineq_rgtside_fnr <- abs(mean(dat$Y*A_col) - mean(dat$Y*A_col_prob))
  
  a1_ineq_leftside_fpr <- abs(mean((1-dat$Y)*dat$Yhat*A_col_prob) - mean((1-dat$Y)*dat$Yhat*A_col))
  a1_ineq_rgtside_fpr <- abs(mean((1-dat$Y)*A_col) - mean((1-dat$Y)*A_col_prob))
  
  a1_ineq_leftside_ppv <- abs(mean(dat$Y*dat$Yhat*A_col_prob) - mean(dat$Y*dat$Yhat*A_col))
  a1_ineq_rgtside_ppv <- abs(mean(dat$Yhat*A_col) - mean(dat$Yhat*A_col_prob))
  
  a1_ineq_leftside_npv <- abs(mean((1-dat$Y)*(1-dat$Yhat)*A_col_prob) - mean((1-dat$Y)*(1-dat$Yhat)*A_col))
  a1_ineq_rgtside_npv <- abs(mean((1-dat$Yhat)*A_col) - mean((1-dat$Yhat)*A_col_prob))
  
  # Sensitivity analysis version of bias
  ## calculate third sensitivity parameter
  sens_param_3_Y1 <- mean(A_col*dat$Y)/mean(dat$Y)
  sens_param_3_Yhat1 <- mean(A_col*dat$Yhat)/mean(dat$Yhat)
  sens_param_3_Y0 <- mean(A_col*(1-dat$Y))/mean(1-dat$Y)
  sens_param_3_Yhat0 <- mean(A_col*(1-dat$Yhat))/mean(1-dat$Yhat)
  ## remove impossible combinations of epsilon, epsilon'
  expec_Y1Yhat1 <- mean(A_col_prob*dat$Y*dat$Yhat)/mean(dat$Y*dat$Yhat)
  expec_Y1Yhat0 <- mean(A_col_prob*dat$Y*(1-dat$Yhat))/mean(dat$Y*(1-dat$Yhat))
  expec_Y0Yhat1 <- mean(A_col_prob*(1-dat$Y)*dat$Yhat)/mean((1-dat$Y)*dat$Yhat)
  expec_Y0Yhat0 <- mean(A_col_prob*(1-dat$Y)*(1-dat$Yhat))/mean((1-dat$Y)*(1-dat$Yhat))
  
  if(epsilon>=expec_Y1Yhat0-1 & epsilon<=expec_Y1Yhat0 & epsilon_prime>=expec_Y1Yhat1-1 & epsilon_prime<=expec_Y1Yhat1) {
    bias_epsilon_fnr <- ((1-fnr_wgt)*fnr_marg*epsilon - fnr_wgt*(1-fnr_marg)*epsilon_prime)/sens_param_3_Y1
    fnr_epsilon <- fnr_wgt+bias_epsilon_fnr
  } else { fnr_epsilon<-NA_real_ }
  
  if(epsilon>=expec_Y1Yhat1-1 & epsilon<=expec_Y1Yhat1 & epsilon_prime>=expec_Y1Yhat0-1 & epsilon_prime<=expec_Y1Yhat0) {
    bias_epsilon_tpr <- ((1-tpr_wgt)*tpr_marg*epsilon - tpr_wgt*(1-tpr_marg)*epsilon_prime)/sens_param_3_Y1
    tpr_epsilon <- tpr_wgt+bias_epsilon_tpr
  } else { tpr_epsilon<-NA_real_ }
  
  if(epsilon>=expec_Y0Yhat1-1 & epsilon<=expec_Y0Yhat1 & epsilon_prime>=expec_Y0Yhat0-1 & epsilon_prime<=expec_Y0Yhat0) {
    bias_epsilon_fpr <- ((1-fpr_wgt)*fpr_marg*epsilon - fpr_wgt*(1-fpr_marg)*epsilon_prime)/sens_param_3_Y0
    fpr_epsilon <- fpr_wgt+bias_epsilon_fpr
  } else { fpr_epsilon<-NA_real_ }
  
  if(epsilon>=expec_Y0Yhat0-1 & epsilon<=expec_Y0Yhat0 & epsilon_prime>=expec_Y0Yhat1-1 & epsilon_prime<=expec_Y0Yhat1) {
    bias_epsilon_tnr <- ((1-tnr_wgt)*tnr_marg*epsilon - tnr_wgt*(1-tnr_marg)*epsilon_prime)/sens_param_3_Y0
    tnr_epsilon <- tnr_wgt+bias_epsilon_tnr
  } else { tnr_epsilon<-NA_real_ }
  
  if(epsilon>=expec_Y1Yhat1-1 & epsilon<=expec_Y1Yhat1 & epsilon_prime>=expec_Y0Yhat1-1 & epsilon_prime<=expec_Y0Yhat1) {
    bias_epsilon_ppv <- ((1-ppv_wgt)*ppv_marg*epsilon - ppv_wgt*(1-ppv_marg)*epsilon_prime)/sens_param_3_Yhat1
    ppv_epsilon <- ppv_wgt+bias_epsilon_ppv
  } else { ppv_epsilon<-NA_real_ }
  
  if(epsilon>=expec_Y0Yhat0-1 & epsilon<=expec_Y0Yhat0 & epsilon_prime>=expec_Y1Yhat0-1 & epsilon_prime<=expec_Y1Yhat0) {
    bias_epsilon_npv <- ((1-npv_wgt)*npv_marg*epsilon - npv_wgt*(1-npv_marg)*epsilon_prime)/sens_param_3_Yhat0
    npv_epsilon <- npv_wgt+bias_epsilon_npv
  } else { npv_epsilon<-NA_real_ }
  
  # Named vector of values to return
  return_vec <- sapply(return_vars_covbias, function(x) { get(x) }, USE.NAMES = T)
  return(return_vec)
}

# Simulate a population and calculate metrics ###########################
## Inputs:
### param_1: conditional dependence parameter
### param_2: mean shift parameter (controls probability ratio)
### param_3: covariance parameter (controls AUC)
### N_pop: population size
### p_Y: Intercept for logit(P(Y=1))

## Outputs:
### data: population data
### pop_results: population metrics
#########################################################################
sim_bias_pop <- function(param_1, param_2, param_3, N_pop, p_Y) {
  ## Simulate Z
  Z <- rnorm(N_pop, -0.4, sd = 1)
  # Simulate A and P(A=a|Z)
  Sigma_mat <- matrix(c(20, param_3, param_3, 20), nrow = 2)
  pA_mat <- t(apply(cbind(Z, Z+param_2), MARGIN = 1, function(r) {
    MASS::mvrnorm(n = 1,
                  mu = r,
                  Sigma = Sigma_mat)
  }))
  A <- rbinom(N_pop, 1, expit(pA_mat[,1]))
  
  dat <- data.frame(
    A = A,
    Pa = expit(pA_mat[, 2]),
    Z = Z
  )
  
  # Simulate X: one component dependent on Z
  dat$X_2 <- rnorm(N_pop, 0, sd = 0.5)
  dat$X_3 <- rnorm(N_pop, 1, sd = 0.5)
  dat$X_4 <- rnorm(N_pop, -1, sd = 0.5)
  
  # Simulate Y
  prob_Y <- expit(p_Y + Z + dat$X_2 + dat$X_3 + dat$X_4 + param_1*A)
  dat$Y <- rbinom(N_pop, 1, prob = prob_Y)
  
  # Split training and test data
  dat_train_ind <- sample(1:nrow(dat), .25*nrow(dat))
  dat_train <- dat[dat_train_ind,]
  dat_test <- dat[-dat_train_ind,]
  
  # Yhat based on model trained on training data
  model_Y <- glm(Y~A + Z + X_3 + X_3 + X_4, data = dat_train, family = "binomial")
  Yhat_prob <- predict(model_Y, newdata = dat_test, type = "response")
  dat_test$Yhat <- if_else(Yhat_prob >= 0.5, 1, 0)
  dat <- dat_test
  N_pop_test <- nrow(dat)
  
  # Get metrics on whole population
  pop_results <- get_avals_bias(dat, epsilon=0.05, epsilon_prime=0.05)
  names(pop_results) <- do.call(paste0, expand.grid(return_vars_covbias, c(".0",".1")))
  
  return(list(data = dat, pop_results = pop_results))
}

# Simulate a population and sample to capture both sampling error and bias #############
## Inputs: 
### param_1: conditional dependence parameter
### param_2: mean shift parameter (controls probability ratio)
### param_3: covariance parameter (controls AUC)
### N_pop: population size
### N_sample: sample size
### p_Y: Intercept for logit(P(Y=1))
### R: number of replications
### epsilon, epsilon': sensitivity parameters for estimating bias

## Outputs:
### if length of epsilon and epsilon' = 1: matrix with R rows and 2*length(return_vars_bias) columns
### if length of either is >1: list with R elements, each element is a matrix with 2*length(return_vars_bias) columns and one row for each combination of epsilon, epsilon'
########################################################################################
sim_bias_sample <- function(data_pop, N_sample, R, epsilon, epsilon_prime) {
  N_pop <- nrow(data_pop)
  
  # Results matrix/list for samples
  if(length(epsilon)==1 & length(epsilon_prime)==1) {
    samp_results <- matrix(nrow = R, ncol = length(return_vars_covbias)*2)
    colnames(samp_results) <- c(do.call(paste0, expand.grid(return_vars_covbias, c(".0",".1"))))
  } else {
    epsilon_grid <- tidyr::expand_grid(epsilon, epsilon_prime)
    samp_results <- vector(mode = "list", length=R)
  }
    
  # Draw samples and calculate metrics
  for(i in 1:R) {
    sample_ind <- sample(1:N_pop, N_sample, replace = F)
    dat_samp <- data_pop[sample_ind,]
    
    ## Get bias terms
    if(length(epsilon)==1 & length(epsilon_prime)==1) {
      bias_res <- get_avals_bias(dat_samp, epsilon=epsilon, epsilon_prime = epsilon_prime)
      samp_results[i,] <- bias_res
    } else {
      e_mat <- matrix(nrow = nrow(epsilon_grid), ncol = length(return_vars_covbias)*2)
      colnames(e_mat) <- c(do.call(paste0, expand.grid(return_vars_covbias, c(".0",".1"))))
      ## Loop over epsilon grid values
      for(j in 1:nrow(epsilon_grid)) {
        e = epsilon_grid[[j,1]]
        e_prime = epsilon_grid[[j,2]]
        bias_res <- get_avals_bias(dat_samp, epsilon=e, epsilon_prime=e_prime)
        e_mat[j,] <- bias_res
      }
      ## Add matrix of epsilon, epsilon' results to replications list
      samp_results[[i]] <- e_mat
    }
  }
  return(samp_results)
}

# Simulate data and get bias terms ####################################################################################
## Inputs:
### param_1; parameter controlling conditional dependence of Y, A given Z (effect of A in P(Y=1))
### param_2: parameter controlling mean of P(A=a|Z), as compared to mean of P(A=a)
### param_3: parameter controlling variance of P(A=a|Z) -- larger value lowers the variance, spreading probabilities more over (0,1)
### N: total sample size
### p_Y: Intercept for logit(P(Y=1))
### R: number of replications
### epsilon, epsilon': sensitivity parameters for estimating bias

## Outputs:
### if length of epsilon and epsilon' = 1: matrix with R rows and 2*length(return_vars_bias) columns
### if length of either is >1: list with R elements, each element is a matrix with 2*length(return_vars_bias) columns and one row for each combination of epsilon, epsilon'
##################################################################################################################
sim_bias <- function(param_1, param_2, param_3, param_4=0, N, p_Y, R,
                             epsilon, epsilon_prime) {
  # Results matrix/list
  if(length(epsilon)==1 & length(epsilon_prime)==1) {
    sim_results <- matrix(nrow = R, ncol = length(return_vars_covbias)*2)
    colnames(sim_results) <- c(do.call(paste0, expand.grid(return_vars_covbias, c(".0",".1"))))
  } else {
    epsilon_grid <- tidyr::expand_grid(epsilon, epsilon_prime)
    sim_results <- vector(mode = "list", length=R)
  }
  # Vector for AUC values
  auc_results <- vector(length = R)
  
  for(i in 1:R) {
    ## Simulate Z
    Z <- rnorm(N, -0.4, sd = 1)
    # Simulate A and P(A=a|Z)
    Sigma_mat <- matrix(c(20, param_3, param_3, 20), nrow = 2)
    pA_mat <- t(apply(cbind(Z, Z+param_2), MARGIN = 1, function(r) {
      MASS::mvrnorm(n = 1,
                    mu = r,
                    Sigma = Sigma_mat)
    }))
    A <- rbinom(N, 1, expit(pA_mat[,1]))
    
    dat <- data.frame(
      A = A,
      Pa = expit(pA_mat[, 2]),
      Z = Z
    )
    test_roc = pROC::roc(response = dat$A, predictor = dat$Pa, quiet = TRUE)
    test_auc = pROC::auc(test_roc)
    
    # Simulate X: one component dependent on Z
    dat$X_2 <- rnorm(N, 0, sd = 0.5)
    dat$X_3 <- rnorm(N, 1, sd = 0.5)
    dat$X_4 <- rnorm(N, -1, sd = 0.5)
    
    # Simulate Y
    prob_Y <- expit(p_Y + Z + dat$X_2 + dat$X_3 + dat$X_4 + param_1*A)
    dat$Y <- rbinom(N, 1, prob = prob_Y)
    
    # Split training and test data
    dat_train_ind <- sample(1:nrow(dat), .25*nrow(dat))
    dat_train <- dat[dat_train_ind,]
    dat_test <- dat[-dat_train_ind,]
    
    # Yhat based on model trained on training data
    model_Y <- glm(Y~A + Z + X_3 + X_3 + X_4, data = dat_train, family = "binomial")
    Yhat_prob <- predict(model_Y, newdata = dat_test, type = "response")
    Yhat_prob <- Yhat_prob + rnorm(nrow(dat_test), mean = 0, sd = param_4)
    dat_test$Yhat <- if_else(Yhat_prob >= 0.5, 1, 0)
    dat <- dat_test
    
    # Get bias terms
    if(length(epsilon)==1 & length(epsilon_prime)==1) {
      bias_res <- get_avals_bias(dat, epsilon=epsilon, epsilon_prime = epsilon_prime)
      sim_results[i,] <- bias_res
    } else {
      e_mat <- matrix(nrow = nrow(epsilon_grid), ncol = length(return_vars_covbias)*2)
      colnames(e_mat) <- c(do.call(paste0, expand.grid(return_vars_covbias, c(".0",".1"))))
      ## Loop over epsilon grid values
      for(j in 1:nrow(epsilon_grid)) {
        e = epsilon_grid[[j,1]]
        e_prime = epsilon_grid[[j,2]]
        bias_res <- get_avals_bias(dat, epsilon=e, epsilon_prime=e_prime)
        e_mat[j,] <- bias_res
      }
      ## Add matrix of epsilon, epsilon' results to replications list
      sim_results[[i]] <- e_mat
    }
    auc_results[i] <- test_auc
  }
  return(list(sim_results = sim_results, auc_results = auc_results))
}

# Shape simulation results to get means, %-tile intervals ########################################################################
## Inputs:
### sim_results: list (R components) or matrix (R rows) of results
### epsilon, epsilon': sensitivity parameters used for estimating bias

## Outputs:
### if length of epsilon and epsilon' = 1: dataframe with columns metric, mean, ci_low (0.025 quantile), ci_high (0.975 quantile)
### if length of either is >1: extra columns for epsilon, epsilon'
##################################################################################################################################
sim_shape <- function(sim_results, epsilon, epsilon_prime) {
  # If sim_results is from sim_bias_sample, extract population metrics
  if("pop_results" %in% names(sim_results)) {
    sim_results_df <- sim_results$samp_results
    pop_vec <- sim_results$pop_results
  } else { sim_results_df <- sim_results$sim_results }
  
  # If single epsilon, epsilon' combination, aggregate over columns. Otherwise, aggregate over matrices.
  if(length(epsilon)==1 & length(epsilon_prime)==1) {
    ## Aggregate over columns: get means, 95%-tile interval
    mean_vec <- apply(sim_results_df, MARGIN=2, mean, na.rm=T)
    cilow_vec <- apply(sim_results_df, MARGIN=2, quantile, probs = 0.025, na.rm=T)
    cihigh_vec <- apply(sim_results_df, MARGIN=2, quantile, probs = 0.975, na.rm=T)
    res <- tibble(metric = names(mean_vec), mean = mean_vec, cilow = cilow_vec, cihigh = cihigh_vec)
    ## Add population metrics if they exist
    if("pop_results" %in% names(sim_results)) {
      res <- left_join(res, tibble(metric = names(pop_vec), pop = pop_vec), by = "metric")
    }
  } else {
    ## Aggregate over matrices: get means, 95%-tile interval
    epsilon_grid <- tidyr::expand_grid(epsilon, epsilon_prime)
    colnames(epsilon_grid) <- c("e_col", "eprime_col")
    mean_mat <- cbind(epsilon_grid, apply(simplify2array(sim_results_df), 1:2, mean, na.rm=T))
    ci_mat <- apply(simplify2array(sim_results_df), 1:2, quantile, prob = c(0.025, 0.975), na.rm=T)
    cilow_mat <- cbind(epsilon_grid, ci_mat[1,,])
    cihigh_mat <- cbind(epsilon_grid, ci_mat[2,,])
    
    mean_long <- pivot_longer(mean_mat, -c(e_col, eprime_col), names_to = "metric", values_to = "mean")
    ## Add population metrics if they exist
    if("pop_results" %in% names(sim_results)) {
      mean_long_pop <- left_join(mean_long, tibble(metric = names(pop_vec), pop = pop_vec), by = "metric")
    } else { mean_long_pop <- mean_long }
    cilow_long <- pivot_longer(cilow_mat, -c(e_col, eprime_col), names_to = "metric", values_to = "cilow")
    cihigh_long <- pivot_longer(cihigh_mat, -c(e_col, eprime_col), names_to = "metric", values_to = "cihigh")
    
    res <- left_join(mean_long_pop, cilow_long, by = c("e_col", "eprime_col", "metric")) %>%
      left_join(cihigh_long, by = c("e_col", "eprime_col", "metric"))
  }
  return(res)
}

# Helper function to aggregate results over a specific subset of epsilon, epsilon' ranges ###############################
## Inputs:
### sim_results_long: simulation results as shaped by sim_shape, with split metric and group columns
### matrics: metric/bias terms to retain in the results (any of return_vars_covbias)
### param_vary: parameter to group results by (e.g., param1, N)
### epsilon_sub: range of epsilon (epsilon') values to retain

## Outputs: summary tibble with columns param_vary, group, and max/min mean, cilow, cihigh columns for each of 'metrics'
#########################################################################################################################
sim_shape_epsilonsub <- function(sim_results_long, metrics, param_vary, epsilon_sub, epsilon_prime_sub) {
  
  sim_results_long %>% mutate(across(c(e_col, eprime_col), \(x) round(x,2))) %>%
    filter(e_col>=min(epsilon_sub), e_col<=max(epsilon_sub), eprime_col>=min(epsilon_prime_sub), eprime_col<=max(epsilon_prime_sub)) %>%
    filter(metric %in% metrics) %>%
    group_by(.data[[param_vary]],metric,group) %>%
    summarize(meanmax = max(mean), meanmin = min(mean),
              cilowmax = max(cilow), cilowmin = min(cilow),
              cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
    pivot_wider(id_cols = c(all_of(param_vary),group), names_from = metric, values_from = c(meanmax,meanmin, cilowmax,cilowmin, cihighmax,cihighmin), names_sep = c("_")) %>%
    mutate(e_range = gsub("0.","", as.character(max(epsilon_sub)), fixed = T))
}
