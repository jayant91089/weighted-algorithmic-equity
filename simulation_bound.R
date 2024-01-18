library(latex2exp)
source('simulation_functions.R')
# Load saved BISG accuracy results
load("rda_files/param_bisg_mean.rda")
load("rda_files/param_bisg_auc.rda")

#####################
# Setup and run sim
#####################
facet_labs_A <- c("0" = "A=0", "1" = "A=1")

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

###### NOTE: Run simulations in simulation_epsilon first, then load saved results here. ####################
load("rda_files/param1_epsilon.rda")
load("rda_files/param2_epsilon.rda")
load("rda_files/param3_epsilon.rda")

# Shape results
param1_res_long <- lapply(param1_res, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param1") %>%
  mutate(param1 = as.numeric(param1),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])
param2_res_long <- lapply(param2_res, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param2") %>%
  mutate(param2 = as.numeric(param2),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])
param3_res_long <- lapply(param3_res, function(r) {
  sim_shape(r, epsilon = epsilon_vals, epsilon_prime = epsilon_prime_vals)
}) %>% bind_rows(.id = "param3") %>%
  mutate(param3 = as.numeric(param3),
         group = str_split(metric, "\\.", simplify = T)[,2],
         metric = str_split(metric, "\\.", simplify = T)[,1])

test_auc_3 <- unlist(lapply(param3_res, function(r) {
  r$auc_results[1]
}))

# Transformed versions for plots
param1_res_plot <- filter(param1_res_long, metric %in% c("bias_bound_abs_fnr", "bias_obs_abs_fnr", "bias_bound_abs_fpr", "bias_obs_abs_fpr",
                                                         "bias_bound_abs_ppv", "bias_obs_abs_ppv", "bias_bound_abs_npv", "bias_obs_abs_npv",
                                                         "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr",
                                                         "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr",
                                                         "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                                                         "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv")) %>%
  group_by(param1,metric,group) %>%
  summarize(meanmax = max(mean), meanmin = min(mean),
            cilowmax = max(cilow), cilowmin = min(cilow),
            cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
  pivot_wider(id_cols = c(param1,group), names_from = metric, values_from = c(meanmax, meanmin, cilowmax, cilowmin, cihighmax, cihighmin), names_sep = c("_"))

param2_res_plot <- filter(param2_res_long, metric %in% c("bias_bound_abs_fnr", "bias_obs_abs_fnr", "bias_bound_abs_fpr", "bias_obs_abs_fpr",
                                                         "bias_bound_abs_ppv", "bias_obs_abs_ppv", "bias_bound_abs_npv", "bias_obs_abs_npv",
                                                         "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr", "pa_ratio",
                                                         "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr",
                                                         "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                                                         "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv")) %>%
  group_by(param2,metric,group) %>%
  summarize(meanmax = max(mean), meanmin = min(mean),
            cilowmax = max(cilow), cilowmin = min(cilow),
            cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
  pivot_wider(id_cols = c(param2,group), names_from = metric, values_from = c(meanmax, meanmin, cilowmax, cilowmin, cihighmax, cihighmin), names_sep = c("_"))

param3_res_plot <- filter(param3_res_long, metric %in% c("bias_bound_abs_fnr", "bias_obs_abs_fnr", "bias_bound_abs_fpr", "bias_obs_abs_fpr",
                                                         "bias_bound_abs_ppv", "bias_obs_abs_ppv", "bias_bound_abs_npv", "bias_obs_abs_npv",
                                                         "a1_ineq_leftside_fnr", "a1_ineq_rgtside_fnr",
                                                         "a1_ineq_leftside_fpr", "a1_ineq_rgtside_fpr",
                                                         "a1_ineq_leftside_ppv", "a1_ineq_rgtside_ppv",
                                                         "a1_ineq_leftside_npv", "a1_ineq_rgtside_npv")) %>%
  group_by(param3,metric,group) %>%
  summarize(meanmax = max(mean), meanmin = min(mean),
            cilowmax = max(cilow), cilowmin = min(cilow),
            cihighmax = max(cihigh), cihighmin = min(cihigh)) %>%
  pivot_wider(id_cols = c(param3,group), names_from = metric, values_from = c(meanmax, meanmin, cilowmax, cilowmin, cihighmax, cihighmin), names_sep = c("_"))

# Combine parameters 2 and 3
param23_list <- list("param2" = rename(param2_res_plot, "param" = "param2"),
                     "param3" = rename(param3_res_plot, "param" = "param3"))
param23_res <- bind_rows(param23_list, .id = "param_set")

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
## Bias bound
p_part1 <- ggplot(filter(param1_res_plot, group=="1"), aes(x = param1)) +
  geom_line(aes(y = meanmax_bias_bound_abs_fnr, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_fnr, linetype = "observed")) +
  labs(linetype = "Bias measure", x = "Conditional dependence parameter", y = "Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

## Assumption
p_part1_a1 <- ggplot(filter(param1_res_plot, group=="1"), aes(x = param1)) +
  geom_line(aes(y = meanmax_a1_ineq_leftside_fnr, linetype = "delta")) +
  geom_line(aes(y = meanmax_a1_ineq_rgtside_fnr, linetype = "deltastar")) +
  labs(x = "Conditional dependence parameter", linetype = "", y = "Mean of 100 replications") +
  scale_linetype_manual(values = c("solid", "dashed"), labels = unname(TeX(c("$|\\hat{\\delta}|$", "$|\\hat{\\delta}^*|$")))) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

## Other metrics
p_part1_fpr <- ggplot(param1_res_plot, aes(x = param1)) +
  geom_line(aes(y = meanmax_bias_bound_abs_fpr, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_fpr, linetype = "observed")) +
  facet_wrap(vars(group), labeller = as_labeller(facet_labs_A)) +
  labs(linetype = "Bias measure", x = "Conditional dependence parameter", y = "FPR Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(0,0.2))

p_part1_ppv <- ggplot(param1_res_plot, aes(x = param1)) +
  geom_line(aes(y = meanmax_bias_bound_abs_ppv, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_ppv, linetype = "observed")) +
  facet_wrap(vars(group), labeller = as_labeller(facet_labs_A)) +
  labs(linetype = "Bias measure", x = "Conditional dependence parameter", y = "PPV Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(0,0.2))

p_part1_npv <- ggplot(param1_res_plot, aes(x = param1)) +
  geom_line(aes(y = meanmax_bias_bound_abs_npv, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_npv, linetype = "observed")) +
  facet_wrap(vars(group), labeller = as_labeller(facet_labs_A)) +
  labs(linetype = "Bias measure", x = "Conditional dependence parameter", y = "NPV Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(0,0.2))

# Part 2: Probability ratio and AUC
## Probability ratio
p_part2_ratio <- ggplot(param2_res_plot, aes(x = param2, y = mean, linetype = group)) +
  geom_line(aes(y = meanmax_pa_ratio)) +
  labs(y = "Ratio", x = unname(TeX("$\\hat{p}_A$ accuracy parameter: mean")), linetype = "Group") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

## AUC
p_part2_auc <- tibble(param3_vals, test_auc_3) %>%
  ggplot(aes(x = param3_vals, y = test_auc_3)) +
  geom_point() +
  labs(y = 'AUC', x = unname(TeX("$\\hat{p}_A$ accuracy parameter: covariance"))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

## Bias bound
p_part2 <- ggplot(filter(param23_res, group == "1"), aes(x = param)) +
  geom_line(aes(y = meanmax_bias_bound_abs_fnr, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_fnr, linetype = "observed")) +
  geom_vline(data = param_bisg, aes(xintercept = param_min), color = "red") +
  geom_vline(data = param_bisg, aes(xintercept = param_max), color = "red") +
  facet_wrap(vars(param_set), scales = "free_x", labeller = as_labeller(c(param2 = "Mean parameter", param3 = "Covariance parameter"))) +
  labs(x = unname(TeX("$\\hat{p}_A$ accuracy parameter")), linetype = "Bias measure", y = "Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

## Assumption
p_part2_a1 <- ggplot(filter(param23_res, group=="1"), aes(x = param)) +
  geom_line(aes(y = meanmax_a1_ineq_leftside_fnr, linetype = "delta")) +
  geom_line(aes(y = meanmax_a1_ineq_rgtside_fnr, linetype = "deltastar")) +
  facet_wrap(vars(param_set), scales = "free_x", labeller = as_labeller(c(param2 = "Mean parameter", param3 = "Covariance parameter"))) +
  labs(x = unname(TeX("$\\hat{p}_A$ accuracy parameter: mean")), linetype = "", y = "Mean of 100 replications") +
  scale_linetype_manual(values = c("solid", "dashed"), labels = unname(TeX(c("$|\\hat{\\delta}|$", "$|\\hat{\\delta}^*|$")))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

## Other metrics
p_part2_fpr <- ggplot(filter(param23_res, group=="1"), aes(x = param)) +
  geom_line(aes(y = meanmax_bias_bound_abs_fpr, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_fpr, linetype = "observed")) +
  facet_wrap(vars(param_set), scales = "free_x", labeller = as_labeller(c(param2 = "Mean parameter", param3 = "Covariance parameter"))) +
  labs(linetype = "Bias measure", x = unname(TeX("$\\hat{p}_A$ accuracy parameter")), y = "FPR Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(0,0.4))

p_part2_ppv <- ggplot(filter(param23_res, group=="1"), aes(x = param)) +
  geom_line(aes(y = meanmax_bias_bound_abs_ppv, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_ppv, linetype = "observed")) +
  facet_wrap(vars(param_set), scales = "free_x", labeller = as_labeller(c(param2 = "Mean parameter", param3 = "Covariance parameter"))) +
  labs(linetype = "Bias measure", x = unname(TeX("$\\hat{p}_A$ accuracy parameter")), y = "PPV Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(0,0.4))

p_part2_npv <- ggplot(filter(param23_res, group=="1"), aes(x = param)) +
  geom_line(aes(y = meanmax_bias_bound_abs_npv, linetype = "bound")) +
  geom_line(aes(y = meanmax_bias_obs_abs_npv, linetype = "observed")) +
  facet_wrap(vars(param_set), scales = "free_x", labeller = as_labeller(c(param2 = "Mean parameter", param3 = "Covariance parameter"))) +
  labs(linetype = "Bias measure", x = unname(TeX("$\\hat{p}_A$ accuracy parameter")), y = "NPV Bias") +
  scale_linetype_manual(values = c("dashed", "solid"), labels = c("Proposed bound", "Observed")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim = c(0,0.4))
