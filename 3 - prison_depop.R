rm(list=ls())
source('0 - prison_fxns.R')

#############################################################
# Impact of decarceration - R proportional to N
#############################################################

R_ARR <- c(1.5,2,3,5,8)
CONTROL_ARR <- seq(0,.5,.01)
PROP_FACTOR_ARR <- seq(0,1)

combo_decarcerate_res <- expand_grid(R = R_ARR, control = CONTROL_ARR, prop_factor = PROP_FACTOR_ARR) %>% group_by(R,control,prop_factor) %>% do({
  
  # Number infectious contacts for prisoners in 100 days
  adj_size <- N_0 * (1-.$control[1])
  avg_days_for_intro <- 1/(PREV*N_C*adj_size*SAR)
  number_inefctious_contacts <- SIM_TIME/avg_days_for_intro
  
  # Prob outbreak in 100 days
  effective_R <- .$R[1] * (1-.$control[1]*.$prop_factor[1])
  prob_ext <- calc_prob_ext(MAX_CONTAINED, effective_R, K_S)
  prob_ob <- 1 - prob_ext^number_inefctious_contacts
  
  # Average number of infected / hospitalized from an outbreak starting within 100 days
  seir_res <- explore_decarceration(prop_capacity_arr = 1-.$control[1], R_max_arr = .$R[1], prop_factor = .$prop_factor[1], unit_size = N_0, time_step = TIME_STEP)
  mean_inf <- prob_ob * seir_res$tot_inf
#  mean_hosp <- prob_ob * seir_res$tot_hosp
  
  tibble(prob_ob = prob_ob, ob_size = seir_res$tot_inf, mean_inf = mean_inf)
})
write.csv(combo_decarcerate_res,"Results/combo_decarcerate_res.csv",row.names = FALSE)
combo_decarcerate_res <- read_csv("Results/combo_decarcerate_res.csv")

#Linear extraploation
simple <- combo_decarcerate_res %>% filter(control == 0) %>% group_by(R) %>% do({
  control_arr <- c(min(CONTROL_ARR), max(CONTROL_ARR))
  simple_inf_arr <- .$mean_inf[1] * (1-control_arr)/ (1-min(CONTROL_ARR))
  tibble(control = control_arr, mean_inf = simple_inf_arr)
})

control_0 <- combo_decarcerate_res %>% filter(control == 0, prop_factor == 1) %>% rename(base_inf = mean_inf, base_prob = prob_ob, base_size = ob_size) %>%
  ungroup() %>% select(R, base_prob, base_size, base_inf)
combo_decarcerate_res <- left_join(combo_decarcerate_res, control_0, by = 'R')
# 
# p1<-ggplot(combo_decarcerate_res %>% filter(prop_factor == 1)) +
#   geom_line(aes(x=control, y = mean_inf, col = as.factor(R))) +
#   geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
#   theme(legend.title = element_text(size = 12),text = element_text(size=20))
# (p1)
#ggsave('Figs/combo_impact.jpg')

p_prob<-ggplot(combo_decarcerate_res) +
  geom_line(aes(x=control, y = prob_ob, col = as.factor(R), linetype = as.factor(prop_factor))) +
  #  geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
  scale_color_discrete(name = 'Baseline R') +
  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
  xlab('Control') + ylab('Probabilty of an an outbreak') +
  guides(linetype = FALSE) + theme(legend.position = 'top') +
  theme(legend.title = element_text(size = 16),text = element_text(size=20))
(p_prob)

p_prob_n<-ggplot(combo_decarcerate_res) +
  geom_line(aes(x=control, y = prob_ob/base_prob, col = as.factor(R),linetype = as.factor(prop_factor))) +
  geom_line(data = simple, aes(x=control, y = 1-control)) +
  scale_color_discrete(name = 'Baseline R') +
  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
  xlab('Control') + ylab('Normalized outbreak probability') +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 1),text = element_text(size=20))
(p_prob_n)

p_size<-ggplot(combo_decarcerate_res) +
  geom_line(aes(x=control, y = ob_size, col = as.factor(R), linetype = as.factor(prop_factor))) +
  #  geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
  scale_color_discrete(name = 'Baseline R') +
  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
  xlab('Control') + ylab('Outbreak size') +
  guides(color = FALSE) + theme(legend.position = 'top') +
  theme(legend.title = element_text(size = 16),text = element_text(size=20))
(p_size)

p_size_n<-ggplot(combo_decarcerate_res) +
  geom_line(aes(x=control, y = ob_size/base_size, col = as.factor(R),linetype = as.factor(prop_factor))) +
  geom_line(data = simple, aes(x=control, y = 1-control)) +
  scale_color_discrete(name = 'Baseline R') +
  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
  xlab('Control') + ylab('Normalized outbreak size') +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 12),text = element_text(size=20))
(p_size_n)

p_tot<-ggplot(combo_decarcerate_res) +
  geom_line(aes(x=control, y = mean_inf, col = as.factor(R), linetype = as.factor(prop_factor))) +
  #  geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
  scale_color_discrete(name = 'Baseline R') +
  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
  xlab('Control') + ylab('Infections expected per 100 days') +
  theme(legend.position = 'top') +
  theme(legend.title = element_text(size = 16),text = element_text(size=20))
(p_tot)

p_tot_n<-ggplot(combo_decarcerate_res) +
  geom_line(aes(x=control, y = mean_inf/base_inf, col = as.factor(R),linetype = as.factor(prop_factor))) +
  geom_line(data = simple, aes(x=control, y = 1-control)) +
  scale_color_discrete(name = 'Baseline R') +
  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
  xlab('Control') + ylab('Normalized infection rate') +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 12),text = element_text(size=20))
(p_tot_n)

#g_impact <- grid.arrange(p_prob,p_size,p_tot,p_prob_n,p_size_n,p_tot_n, nrow = 2)
g_impact <- grid.arrange(p_tot,p_tot_n, nrow = 2)
ggsave('Figs/depop_impact.jpg',plot = g_impact, scale = 1.2)
g_component <- grid.arrange(p_prob,p_size,p_prob_n,p_size_n, nrow = 2)
ggsave('Figs/depop_component.jpg',plot = g_component, scale = 1.25)

# Same sample numbers
control_0 <- combo_decarcerate_res %>% filter (control == 0, prop_factor == 1)
control_20 <- combo_decarcerate_res %>% filter (control == 0.2, prop_factor == 1)
1- control_20$mean_inf/control_0$mean_inf
