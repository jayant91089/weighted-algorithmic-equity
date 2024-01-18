library(latex2exp)
source('simulation_functions.R')
# Load saved BISG accuracy results
load("rda_files/param_bisg_mean.rda")
load("rda_files/param_bisg_auc.rda")

#####################
# Setup and run sim
#####################
param1_vals <- seq(-0.5, 0.5, 0.1)
param2_vals <- seq(-1, 1, 0.1)
param3_vals <- seq(0, 20, 1)

names(param1_vals) <- param1_vals
names(param2_vals) <- param2_vals
names(param3_vals) <- param3_vals

N <- 5000
p_Y <- c(-0.2)

epsilon_vals <- c(seq(-0.05,0.05,0.01))
epsilon_prime_vals <- c(seq(-0.05,0.05,0.01))

# Parameter 1
param1_res <- lapply(param1_vals, function(p) {
  sim_bias(param_1=p, param_2=0, param_3=20, N=N, p_Y=p_Y, R = 100,
           epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
})
save(param1_res, file = "rda_files/param1_epsilon.rda")

param1_res_long <- lapply(param1_res, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param1") %>%
  mutate(param1 = as.numeric(param1),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])

## Subsets of epsilon_grid
param1_res_sub_02 <- sim_shape_epsilonsub(param1_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param1",
                                          epsilon_sub=c(-0.02,0.02), epsilon_prime_sub=c(-0.02,0.02))
param1_res_sub_03 <- sim_shape_epsilonsub(param1_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param1",
                                          epsilon_sub=c(-0.03,0.03), epsilon_prime_sub=c(-0.03,0.03))
param1_res_sub_04 <- sim_shape_epsilonsub(param1_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param1",
                                          epsilon_sub=c(-0.04,0.04), epsilon_prime_sub=c(-0.04,0.04))

## Transformed version for plot
param1_res_plot <- filter(param1_res_long, metric %in% c("bias_obs_fnr", "bias_epsilon_fnr")) %>%
  group_by(param1,metric,group) %>%
  summarize(meanmax = max(mean), meanmin = min(mean),
            cilowmax = max(cilow), cilowmin = min(cilow),
            cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
  pivot_wider(id_cols = c(param1,group), names_from = metric, values_from = c(meanmax, meanmin, cilowmax, cilowmin, cihighmax, cihighmin), names_sep = c("_")) %>%
  mutate(e_range = "05") %>%
  bind_rows(param1_res_sub_02) %>%
  bind_rows(param1_res_sub_03) %>%
  bind_rows(param1_res_sub_04) %>%
  ### manually set alpha values for epsilon, epsilon' ranges
  mutate(e_range = case_when(
    e_range == "02" ~ 0.6,
    e_range == "03" ~ 0.5,
    e_range == "04" ~ 0.3,
    e_range == "05" ~ 0.1
  ))

# Parameter 2
param2_res <- lapply(param2_vals, function(p) {
  sim_bias(param_1 = 0.25, param_2=p, param_3=20, N=N, p_Y=p_Y, R = 100,
                 epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
})
save(param2_res, file = "rda_files/param2_epsilon.rda")

param2_res_long <- lapply(param2_res, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param2") %>%
  mutate(param2 = as.numeric(param2),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])

## Subsets of epsilon_grid
param2_res_sub_02 <- sim_shape_epsilonsub(param2_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param2",
                                          epsilon_sub=c(-0.02,0.02), epsilon_prime_sub=c(-0.02,0.02))
param2_res_sub_03 <- sim_shape_epsilonsub(param2_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param2",
                                          epsilon_sub=c(-0.03,0.03), epsilon_prime_sub=c(-0.03,0.03))
param2_res_sub_04 <- sim_shape_epsilonsub(param2_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param2",
                                          epsilon_sub=c(-0.04,0.04), epsilon_prime_sub=c(-0.04,0.04))

## Transformed version for plot
param2_res_plot <- filter(param2_res_long, metric %in% c("bias_obs_fnr", "bias_epsilon_fnr")) %>%
  group_by(param2,metric,group) %>%
  summarize(meanmax = max(mean), meanmin = min(mean),
            cilowmax = max(cilow), cilowmin = min(cilow),
            cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
  pivot_wider(id_cols = c(param2,group), names_from = metric, values_from = c(meanmax, meanmin, cilowmax, cilowmin, cihighmax, cihighmin), names_sep = c("_")) %>%
  mutate(e_range = "05") %>%
  bind_rows(param2_res_sub_02) %>%
  bind_rows(param2_res_sub_03) %>%
  bind_rows(param2_res_sub_04) %>%
  ### manually set alpha values for epsilon, epsilon' ranges
  mutate(e_range = case_when(
    e_range == "02" ~ 0.6,
    e_range == "03" ~ 0.5,
    e_range == "04" ~ 0.3,
    e_range == "05" ~ 0.1
  ))

# Parameter 3
param3_res <- lapply(seq_along(param3_vals), function(i) {
  sim_bias(param_1=0.25, param_2=0, param_3 = param3_vals[i], N=N, p_Y=p_Y, R = 100,
            epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
})
save(param3_res, file = "rda_files/param3_epsilon.rda")

param3_res_long <- lapply(param3_res, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param3") %>%
  mutate(param3 = as.numeric(param3),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])

test_auc_3 <- unlist(lapply(param3_res, function(r) {
  r$auc_results[1]
}))

## Subsets of epsilon_grid
param3_res_sub_02 <- sim_shape_epsilonsub(param3_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param3",
                                          epsilon_sub=c(-0.02,0.02), epsilon_prime_sub=c(-0.02,0.02))
param3_res_sub_03 <- sim_shape_epsilonsub(param3_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param3",
                                          epsilon_sub=c(-0.03,0.03), epsilon_prime_sub=c(-0.03,0.03))
param3_res_sub_04 <- sim_shape_epsilonsub(param3_res_long, metrics=c("bias_obs_fnr", "bias_epsilon_fnr"), param_vary="param3",
                                          epsilon_sub=c(-0.04,0.04), epsilon_prime_sub=c(-0.04,0.04))

## Transformed version for plot
param3_res_plot <- filter(param3_res_long, metric %in% c("bias_obs_fnr", "bias_epsilon_fnr")) %>%
  group_by(param3,metric,group) %>%
  summarize(meanmax = max(mean), meanmin = min(mean),
            cilowmax = max(cilow), cilowmin = min(cilow),
            cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
  pivot_wider(id_cols = c(param3,group), names_from = metric, values_from = c(meanmax, meanmin, cilowmax, cilowmin, cihighmax, cihighmin), names_sep = c("_")) %>%
  mutate(e_range = "05") %>%
  bind_rows(param3_res_sub_02) %>%
  bind_rows(param3_res_sub_03) %>%
  bind_rows(param3_res_sub_04) %>%
  ### manually set alpha values for epsilon, epsilon' ranges
  mutate(e_range = case_when(
    e_range == "02" ~ 0.6,
    e_range == "03" ~ 0.5,
    e_range == "04" ~ 0.3,
    e_range == "05" ~ 0.1
  ))

# Look up BISG AUC values and get corresponding parameter 3
param3_auc_df <- tibble(param3 = param3_vals, auc = test_auc_3)
param_cov <- param3_auc_df$param3[findInterval(param_auc, param3_auc_df$auc)]

param_bisg <- data.frame(
  param_set = c("param2", "param3"),
  param_min = c(min(param_mean), min(param_cov)),
  param_max = c(max(param_mean), max(param_cov))
)

########
# Plots
########
# Part 1: Conditional dependence
p_param1_epsilon <- ggplot(filter(param1_res_plot, group=="1"),
                           aes(x = param1, group = e_range)) +
  geom_ribbon(aes(ymin = meanmin_bias_epsilon_fnr, ymax = meanmax_bias_epsilon_fnr, alpha = e_range)) +
  geom_line(aes(y = meanmax_bias_obs_fnr, linetype = "observed")) +
  scale_linetype_manual(values = c("solid"), labels = c("Observed bias")) +
  scale_alpha(breaks = c(0.6,0.5,0.3,0.1), range=c(0.1,0.6), labels=c(unname(TeX("$\\{-0.02,0.02\\}$")), unname(TeX("$\\{-0.03,0.03\\}$")), unname(TeX("$\\{-0.04,0.04\\}$")), unname(TeX("$\\{-0.05,0.05\\}$")))) +
  labs(linetype = "", alpha = unname(TeX("$\\epsilon,\\epsilon'$ range")), x = "Conditional dependence parameter", y = "FNR") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(-0.06,0.06))


# Part 2-3: Probability ratio and AUC
p_param2_epsilon <- ggplot(filter(param2_res_plot, group=="1"),
                           aes(x = param2, group = e_range)) +
  geom_ribbon(aes(ymin = meanmin_bias_epsilon_fnr, ymax = meanmax_bias_epsilon_fnr, alpha = e_range)) +
  geom_line(aes(y = meanmax_bias_obs_fnr, linetype = "observed")) +
  geom_vline(data = filter(param_bisg, param_set=="param2"), aes(xintercept = param_min), color = "red") +
  geom_vline(data = filter(param_bisg, param_set=="param2"), aes(xintercept = param_max), color = "red") +
  scale_linetype_manual(values = c("solid"), labels = c("Observed bias")) +
  scale_alpha(breaks = c(0.6,0.5,0.3,0.1), range=c(0.1,0.7), labels=c(unname(TeX("$\\{-0.02,0.02\\}$")), unname(TeX("$\\{-0.03,0.03\\}$")), unname(TeX("$\\{-0.04,0.04\\}$")), unname(TeX("$\\{-0.05,0.05\\}$")))) +
  labs(linetype = "", alpha = unname(TeX("$\\epsilon,\\epsilon'$ range")), x = unname(TeX("$\\hat{p}_A$ accuracy parameter: mean")), y = "FNR") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(-0.06,0.06))

p_param3_epsilon <- ggplot(filter(param3_res_plot, group=="1"),
                           aes(x = param3, group = e_range)) +
  geom_ribbon(aes(ymin = meanmin_bias_epsilon_fnr, ymax = meanmax_bias_epsilon_fnr, alpha = e_range)) +
  geom_line(aes(y = meanmax_bias_obs_fnr, linetype = "observed")) +
  geom_vline(data = filter(param_bisg, param_set=="param3"), aes(xintercept = param_min), color = "red") +
  geom_vline(data = filter(param_bisg, param_set=="param3"), aes(xintercept = param_max), color = "red") +
  scale_linetype_manual(values = c("solid"), labels = c("Observed bias")) +
  scale_alpha(breaks = c(0.6,0.5,0.3,0.1), range=c(0.1,0.7), labels=c(unname(TeX("$\\{-0.02,0.02\\}$")), unname(TeX("$\\{-0.03,0.03\\}$")), unname(TeX("$\\{-0.04,0.04\\}$")), unname(TeX("$\\{-0.05,0.05\\}$")))) +
  labs(linetype = "", alpha = unname(TeX("$\\epsilon,\\epsilon'$ range")), x = unname(TeX("$\\hat{p}_A$ accuracy parameter: covariance")), y = "FNR") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(-0.06,0.06))

## Combine parts 2-3
param23_list <- list("param2" = rename(param2_res_plot, "param" = "param2"),
                     "param3" = rename(param3_res_plot, "param" = "param3"))
param23_res <- bind_rows(param23_list, .id = "param_set")

p_param23_epsilon <- ggplot(filter(param23_res, group=="1"),
                            aes(x = param, group = e_range)) +
  geom_ribbon(aes(ymin = meanmin_bias_epsilon_fnr, ymax = meanmax_bias_epsilon_fnr, alpha = e_range)) +
  geom_line(aes(y = meanmax_bias_obs_fnr, linetype = "observed")) +
  geom_vline(data = param_bisg, aes(xintercept = param_min), color = "red") +
  geom_vline(data = param_bisg, aes(xintercept = param_max), color = "red") +
  facet_wrap(vars(param_set), scales = "free_x", labeller = as_labeller(c(param2 = "Mean parameter", param3 = "Covariance parameter"))) +
  labs(x = unname(TeX("$\\hat{p}_A$ accuracy parameter")), linetype = "", alpha = unname(TeX("$\\epsilon,\\epsilon'$ range")),
       y = "FNR") +
  scale_linetype_manual(values = c("solid"), labels = c("Observed bias")) +
  scale_alpha(breaks = c(0.6,0.5,0.3,0.1), range=c(0.1,0.7), labels=c(unname(TeX("$\\{-0.02,0.02\\}$")), unname(TeX("$\\{-0.03,0.03\\}$")), unname(TeX("$\\{-0.04,0.04\\}$")), unname(TeX("$\\{-0.05,0.05\\}$")))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(-0.06,0.06))
