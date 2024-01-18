library(tidyverse)

########
# Setup
########
fnr_wgt <- c(0.35, 0.48)
mean_fn_overall <- 0.16
mean_pos_overall <- 0.42
fnr_overall <- mean_fn_overall/mean_pos_overall
mean_tp_overall <- mean_pos_overall-mean_fn_overall
mean_est_fn <- c(0.09, 0.07)
mean_est_pos <- mean_est_fn / fnr_wgt
mean_est_tp <- mean_est_pos-mean_est_fn
prop_A <- c(0.48, 0.52)

# Unknown parameter grid
mean_true_fn_vec <- seq(0.005, 1, 0.005)
mean_true_pos_vec <- seq(0.005, 1, 0.005)
mean_true_grid <- expand_grid(mean_true_fn_vec, mean_true_pos_vec)
mean_true_grid <- filter(mean_true_grid, mean_true_fn_vec <= mean_true_pos_vec)
  
bias_mat <- matrix(nrow = nrow(mean_true_grid), ncol = 2)

group_colors <- c("#440154", "#2A788E", "#22A884", "#7AD151")
group_labels <- c("1" = "Group 1", "2" = "Group 2")

#################################################
# Loop over groups and unknown parameter values
#################################################
for(i in 1:2) {
  for(j in 1:nrow(mean_true_grid)) {
    ## Set unknown parameters
    mean_true_fn <- mean_true_grid[j,1]
    mean_true_pos <- mean_true_grid[j,2]
    ## Calculate bias
    bias_val <- (mean_est_fn[i] - mean_true_fn)/mean_true_pos + fnr_wgt[i]*(1-mean_est_pos[i]/mean_true_pos)
    
    bias_mat[j,i] <- bias_val[1,1]
  }
}
bias_mat_red <- as.data.frame(bias_mat)
colnames(bias_mat_red) <- c("G_1", "G_2")

# Remove impossible values (higher than overall group proportions OR mean_true_fn > mean_true_pos)
for(i in 1:2) {
  temp_df <- mean_true_grid %>% mutate(is_impossible = mean_true_fn_vec > mean_fn_overall | mean_true_pos_vec > mean_pos_overall | mean_true_fn_vec > mean_true_pos_vec,
                                       epsilon = (mean_est_fn[i] - mean_true_fn_vec)/mean_fn_overall,
                                       epsilon_prime = (mean_est_tp[i] - (mean_true_pos_vec-mean_true_fn_vec))/mean_tp_overall,
                                       delta = mean_est_fn[i] - mean_true_fn_vec,
                                       delta_star = mean_est_pos[i] - mean_true_pos_vec,
                                       a1_epsilon_rgtside = abs(epsilon)*2*mean_fn_overall/mean_tp_overall,
                                       is_out_epsilon = (mean_est_fn[i]/mean_fn_overall - 1 > epsilon) | (mean_est_fn[i]/mean_fn_overall < epsilon),
                                       is_out_epsilon_prime = (mean_est_tp[i]/mean_tp_overall - 1 > epsilon_prime) | (mean_est_tp[i]/mean_tp_overall < epsilon_prime),
                                       sens_param_3 = mean_true_pos_vec/mean_pos_overall
  )
  ## Remove impossible parameter combinations
  bias_mat_red[,i] <- if_else(temp_df$is_impossible, NA_real_, bias_mat_red[,i])
  ## Mark parameter values outside +/- delta range
  bias_mat_red[[paste0('epsilonout_',i)]] <- if_else(temp_df$is_out_epsilon, TRUE, FALSE)
  bias_mat_red[[paste0('epsilonpout_',i)]] <- if_else(temp_df$is_out_epsilon_prime, TRUE, FALSE)
  ## Columns for epsilon, epsilon'
  bias_mat_red[[paste0('epsilon_',i)]] <- temp_df$epsilon
  bias_mat_red[[paste0('epsilonp_',i)]] <- temp_df$epsilon_prime
  bias_mat_red[[paste0('delta_',i)]] <- temp_df$delta
  bias_mat_red[[paste0('deltastar_',i)]] <- temp_df$delta_star
  bias_mat_red[[paste0('a1rgtside_',i)]] <- temp_df$a1_epsilon_rgtside
  bias_mat_red[[paste0('biasep_',i)]] <- ((1-fnr_wgt[i])*fnr_overall*temp_df$epsilon - fnr_wgt[i]*(1-fnr_overall)*temp_df$epsilon_prime)/temp_df$sens_param_3
  bias_mat_red[[paste0('bias0slope_',i)]] <- ((1-fnr_wgt[i])*fnr_overall)/(fnr_wgt[i]*(1-fnr_overall))
}
bias_df <- cbind(mean_true_grid, bias_mat_red)

# Long df
bias_df_long <- bias_df %>% pivot_longer(-c(mean_true_fn_vec, mean_true_pos_vec), names_to = c("metric", "G"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = metric, values_from = value, names_repair = "unique") %>%
  rename(G = G...3, bias = G...4, bias_ep = biasep, out_epsilon = epsilonout, out_epsilon_prime = epsilonpout, 
         epsilon_prime = epsilonp, delta_star = deltastar, a1_epsilon_rgtside = a1rgtside) %>%
  filter(!is.na(bias))

bias_df_long_minmax <- filter(bias_df_long, out_epsilon==0, out_epsilon_prime==0) %>%
  group_by(mean_true_pos_vec, G) %>%
  mutate(
    min_fn_vary = min(bias),
    max_fn_vary = max(bias)
  ) %>% ungroup() %>%
  group_by(mean_true_fn_vec, G) %>%
  mutate(
    min_pos_vary = min(bias),
    max_pos_vary = max(bias)
  ) %>% ungroup() 

bias_df_long_epsilon <- filter(bias_df_long, out_epsilon==0, out_epsilon_prime==0) %>%
  mutate(a1_met = abs(delta) <= abs(delta_star))

########
# Plots
########
sc_bias_invert <- colorspace::scale_colour_continuous_divergingx(palette = "Cividis", limits=c(min(fnr_wgt)-1,max(fnr_wgt)))

# Contour plot: color denotes bias
p_sens_demo_biascolor <- ggplot(bias_df_long_epsilon, aes(x = round(epsilon,4), y = round(epsilon_prime,4), z = bias, group = G, color = after_stat(level))) +
  geom_contour() +
  geom_abline(aes(slope = bias0slope, intercept = 0), color = "red") +
  facet_wrap(vars(G), nrow = 2, labeller = as_labeller(group_labels)) +
  labs(y = unname(latex2exp::TeX("$\\epsilon'$")), 
       x = unname(latex2exp::TeX("$\\epsilon$")), color = "Estimated Bias") +
  metR::geom_text_contour(mapping = aes(label = round(after_stat(level),2)), stroke = 0.3,
                          size = 3, min.size = 1, skip = 1, label.placer = metR::label_placer_fraction(frac = 0.1)) +
  sc_bias_invert +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Contour plot with region violating assumption 1 shaded
bias_df_G1 <- filter(bias_df_long_epsilon, G=="1")
bias_df_G2 <- filter(bias_df_long_epsilon, G=="2")

dat_poly <- data.frame(G = c(rep("1",3), rep("2",3), rep("1",3), rep("2",3)),
                           id = c(rep(1,6), rep(2,6)),
                       x = c(min(bias_df_G1$epsilon), min(bias_df_G1$epsilon), 0,
                             min(bias_df_G2$epsilon), min(bias_df_G2$epsilon), 0,
                             max(bias_df_G1$epsilon), max(bias_df_G1$epsilon), 0,
                             max(bias_df_G2$epsilon), max(bias_df_G2$epsilon), 0),
                       y = c(max(bias_df_G1$a1_epsilon_rgtside), 0,0,
                             max(bias_df_G2$a1_epsilon_rgtside), 0,0,
                             0, -max(bias_df_G1$a1_epsilon_rgtside), 0,
                             0, -max(bias_df_G2$a1_epsilon_rgtside), 0))

p_sens_demo_shaded <- ggplot(bias_df_long_epsilon, aes(x = round(epsilon,4), y = round(epsilon_prime,4))) +
  geom_contour(aes(color = after_stat(level), z = bias, group = G)) +
  geom_polygon(data = dat_poly, mapping = aes(x=x, y =y, group=id), fill = "gray", alpha = 0.8) +
  facet_wrap(vars(G), nrow = 2, labeller = as_labeller(group_labels)) +
  labs(y = unname(latex2exp::TeX("$\\epsilon'$")), 
       x = unname(latex2exp::TeX("$\\epsilon$")), color = "Estimated Bias") +
  metR::geom_text_contour(mapping = aes(label = round(ggplot2::after_stat(level),2), color = after_stat(level), z = bias, group = G), 
                          stroke = 0.3, size = 3, min.size = 1, skip = 1, label.placer = metR::label_placer_fraction(frac = 0.1)) +
  sc_bias_invert +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(xlim = c(-0.6,0.55), ylim = c(-0.7,max(bias_df_long_epsilon$epsilon_prime)))
