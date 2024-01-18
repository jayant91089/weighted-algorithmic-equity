library(latex2exp)
source('simulation_functions.R')

#####################
# Setup and run sim
#####################
N_vals <- seq(1000,10000, 1000)
param1_vals <- seq(-0.5, 0.5, 0.1)
names(N_vals) <- N_vals
names(param1_vals) <- param1_vals

p_Y <- c(-0.2)

epsilon_vals <- c(seq(-0.05,0.05,0.01))
epsilon_prime_vals <- c(seq(-0.05,0.05,0.01))

################ Sample size
simN_res_pop <- sim_bias_pop(param_1=0.25, param_2=0, param_3=20, N_pop = 50000, p_Y=p_Y)

simN_res_all <- lapply(N_vals, function(n) {
  sim_bias_sample(data_pop = simN_res_pop$data, N_sample=n, R = 100,
                  epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals) 
})

simN_res_all <- lapply(simN_res_all, function(r) {
  list(samp_results = r, pop_results = simN_res_pop$pop_results)
})

simN_res_long <- lapply(simN_res_all, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "N") %>%
  mutate(N = as.numeric(N),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])

## Transformed version for plot
simN_bias_plot <- filter(simN_res_long, metric %in% c("bias_obs_fnr"), e_col==0.05,eprime_col==0.05) %>%
  select(-c(pop, e_col,eprime_col))

################ Conditional dependence: N_sample = 2,000
simparam1_res_all_2K <- lapply(param1_vals, function(p) {
  res_pop <- sim_bias_pop(param_1=p, param_2 = 0, param_3=20, N_pop = 50000, p_Y=p_Y)
  
  res_samp <- sim_bias_sample(data_pop = res_pop$data, N_sample=2000, R = 100, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
  list(samp_results = res_samp, pop_results = res_pop$pop_results)
})

simparam1_res_long_2K <- lapply(simparam1_res_all_2K, function(r) {
  res <- sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param1") %>%
  mutate(param1 = as.numeric(param1),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])

## Transformed version for plot
simparam1_bias_plot_2K <- filter(simparam1_res_long_2K, metric %in% c("bias_obs_fnr"), e_col==0.05,eprime_col==0.05) %>%
  select(-c(pop, e_col,eprime_col))

#########
# Plots
#########
# Sample size
p_N_bias <- ggplot(filter(simN_bias_plot, group=="1"), aes(x = N, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), width = 500) +
  geom_hline(yintercept=0, color = "black", linetype = "solid") +
  labs(x = "Sample size", y = "Bias") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Conditional dependence
p_param1_bias <- ggplot(filter(simparam1_bias_plot_2K, group=="1"), aes(x = param1, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), width = 0.05) +
  geom_hline(yintercept=0, color = "black", linetype = "solid") +
  labs(x = "Conditional dependence parameter", y = "Bias") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
