library(here)
rm(list=ls())
source(here('0 - prison_fxns.R'))

#############################################################
# Impact of decarceration - R proportional to N
#############################################################

R_ARR <- c(1.5,2,3,5,8)
CONTROL_ARR <- seq(0,.5,.01)
MECHANISM_TABLE <- tribble(
  ~mech_num, ~"Mechanism of control", ~N_factor, ~R_factor, ~rho_factor,
  1, "Depopulation\n(freq-dep)", 1, 0,0,
  2, "Vaccination\n(freq-dep)", 0, 0, 1,
  3, "Density-dep\ncontrol", 1, 1, 0,
  4, "Reduced\ntransmission only", 0, 1, 0 
)

combo_decarcerate_res <- expand_grid(R = R_ARR, control = CONTROL_ARR, MECHANISM_TABLE) %>% group_by(R,control,mech_num) %>% do({
  # Number infectious contacts for prisoners in 100 days
  adj_size <- N_0 * (1-.$control*.$N_factor)
  avg_days_for_intro <- 1/(PREV*N_C*adj_size*SAR)
  number_inefctious_contacts <- SIM_TIME/avg_days_for_intro
  
  # Prob outbreak in 100 days
  effective_R <- .$R * (1-.$control*.$R_factor[1])
  prob_ext <- calc_prob_ext(MAX_CONTAINED, effective_R, K_S)
  prob_ob <- 1 - prob_ext^number_inefctious_contacts
  
  # Average number of infected / hospitalized from an outbreak starting within 100 days
  seir_res <- explore_decarceration(prop_capacity_arr = 1, R_max_arr = effective_R, prop_factor = 1, unit_size = adj_size, time_step = TIME_STEP)
  mean_inf <- prob_ob * seir_res$tot_inf * PROP_SYMPTOMATIC * (1-.$control * .$rho_factor)

  tibble(prob_ob = prob_ob, ob_size = seir_res$tot_inf, mean_inf = mean_inf)
})
combo_decarcerate_res <- left_join(combo_decarcerate_res,MECHANISM_TABLE)
#write.csv(combo_decarcerate_res,"Results/combo_decarcerate_res.csv",row.names = FALSE)
combo_decarcerate_res <- read_csv("Results/combo_decarcerate_res.csv")

#Linear extraploation
simple <- combo_decarcerate_res %>% filter(control == 0) %>% group_by(R) %>% do({
  control_arr <- c(min(CONTROL_ARR), max(CONTROL_ARR))
  simple_inf_arr <- .$mean_inf[1] * (1-control_arr)/ (1-min(CONTROL_ARR))
  tibble(control = control_arr, mean_inf = simple_inf_arr)
})

control_0 <- combo_decarcerate_res %>% filter(control == 0) %>% group_by(R) %>% summarize(base_inf = mean(mean_inf), base_prob = mean(prob_ob), base_size = mean(ob_size))
combo_decarcerate_res <- left_join(combo_decarcerate_res, control_0, by = 'R')
# 
# p1<-ggplot(combo_decarcerate_res %>% filter(prop_factor == 1)) +
#   geom_line(aes(x=control, y = mean_inf, col = as.factor(R))) +
#   geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
#   theme(legend.title = element_text(size = 12),text = element_text(size=20))
# (p1)
#ggsave('Figs/combo_impact.jpg')

sparse_control <- CONTROL_ARR[seq(2,length(CONTROL_ARR),5)]
combo_decarcerate_res_sparse <- combo_decarcerate_res %>% filter(control %in% sparse_control)

man_color <- gg_color_hue(length(R_ARR))
names(man_color) <- R_ARR
man_shape <- c(16,2, 8, 6, 15, 0:25)
man_shape <- man_shape[1:length(R_ARR)]
names(man_shape) <- R_ARR
man_line <- MECHANISM_TABLE$mech_num[c(2, 1, 3, 4)]
names(man_line) <- MECHANISM_TABLE$mech_num

###############################

p_8 <- ggplot(combo_decarcerate_res %>% filter (R == 8)) +
  geom_line(aes(x=control, y = mean_inf, col = as.factor(R), linetype = as.factor(mech_num)), size = 1) +
  # geom_point(data = combo_decarcerate_res_sparse %>% filter (R == 8), aes(x=control, y = mean_inf, col =  as.factor(R), shape = as.factor(R), linetype = as.factor(mech_num)), size = 2.5) +
  # geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Expected cases\nper 100 days') + labs(tag = "A)") + expand_limits(y=0) +
  guides(color = "none", shape = "none") + theme(legend.position = "top") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_8)

p_1p5 <- ggplot(combo_decarcerate_res %>% filter (R == 1.5)) +
  geom_line(aes(x=control, y = mean_inf, col = as.factor(R), linetype = as.factor(mech_num)), size = 1) +
  # geom_point(data = combo_decarcerate_res_sparse %>% filter (R == 1.5), aes(x=control, y = mean_inf, col = as.factor(R), shape = as.factor(R), linetype = as.factor(mech_num)), size = 2.5) +
  # geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Expected cases\nper 100 days') + labs(tag = "B)") + expand_limits(y=0) +
  guides(color = "none", shape = "none") + theme(legend.position = "top") +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_1p5)

p_control_n<-
  ggplot(combo_decarcerate_res %>% filter (mech_num %in% c(1,3))) +
  geom_line(data = simple, aes(x=control, y = 1-control), size = 1) +
  geom_line(aes(x=control, y = mean_inf/base_inf, col = as.factor(R),linetype = as.factor(mech_num)), size = 1) +
  geom_point(data = combo_decarcerate_res_sparse %>% filter (mech_num %in% c(3)), aes(x=control, y = mean_inf/base_inf, col = as.factor(R), shape = as.factor(R),linetype = as.factor(mech_num)), size = 2.5) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized\ncase rate') + labs(tag = "C)") +
  guides(linetype = "none") + theme(legend.position = 'top') +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_control_n)

p_R_n<-
  ggplot(combo_decarcerate_res %>% filter (mech_num %in% c(4))) +
  geom_line(data = simple, aes(x=control, y = 1-control), size = 1) +
  geom_line(aes(x=control, y = mean_inf/base_inf, col = as.factor(R),linetype = as.factor(mech_num)), size = 1) +
  geom_point(data = combo_decarcerate_res_sparse %>% filter (mech_num %in% c(4)), aes(x=control, y = mean_inf/base_inf, col = as.factor(R), shape = as.factor(R),linetype = as.factor(mech_num)), size = 2.5) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized\ncase rate') + labs(tag = "D)") +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_R_n)

g_impact <- grid.arrange(p_8, p_control_n, p_1p5, p_R_n, nrow = 2)
ggsave('Figs/depop_impact.jpg',plot = g_impact, scale = 1.25)
ggsave("Figs/TIFF/Figure 2.eps",plot = g_impact, width=4320, units="px", dpi=300)
###############################

p_prob_ctrl<-
  ggplot(combo_decarcerate_res %>% filter (mech_num %in% c(1,3))) +
  geom_line(data = simple, aes(x=control, y = 1-control), color = 'grey', size = 1) +
  geom_line(data = simple, aes(x=control, y = 1), size = 1) +
  geom_line(aes(x=control, y = prob_ob/base_prob, col = as.factor(R),linetype = as.factor(mech_num)), size = 1) +
  geom_point(data = combo_decarcerate_res_sparse %>% filter (mech_num %in% c(3)), aes(x=control, y = prob_ob/base_prob, col = as.factor(R),  shape = as.factor(R),linetype = as.factor(mech_num)), size = 2.5) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized outbreak\nprobability') + labs(tag = "A)") + expand_limits(y=0) +
  guides(color = "none", shape = "none") + theme(legend.position = "top") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_prob_ctrl)

p_prob_R<-
  ggplot(combo_decarcerate_res %>% filter (mech_num %in% c(4))) +
  geom_line(data = simple, aes(x=control, y = 1-control), color = 'grey', size = 1) +
  geom_line(aes(x=control, y = prob_ob/base_prob, col = as.factor(R),linetype = as.factor(mech_num)), size = 1) +
  geom_point(data = combo_decarcerate_res_sparse %>% filter (mech_num %in% c(4)), aes(x=control, y = prob_ob/base_prob, col = as.factor(R),  shape = as.factor(R),linetype = as.factor(mech_num)), size = 2.5) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized outbreak\nprobability') + labs(tag = "B)") +  expand_limits(y=0) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_prob_R)

p_size_ctrl<-ggplot(combo_decarcerate_res %>% filter (mech_num %in% c(3))) +
  geom_line(aes(x=control, y = mean_inf/prob_ob/base_size/PROP_SYMPTOMATIC, col = as.factor(R),linetype = as.factor(mech_num)), size = 1) +
  geom_line(data = simple, aes(x=control, y = 1-control), size = 1) +
  # geom_line(data = simple, aes(x=control, y = .9925-control), linetype = 'dashed', size = 1) +
  geom_point(data = combo_decarcerate_res_sparse %>% filter (mech_num %in% c(3)), aes(x=control, y = mean_inf/prob_ob/base_size/PROP_SYMPTOMATIC, col = as.factor(R), shape = as.factor(R),linetype = as.factor(mech_num)), size = 2.5) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized outbreak size') + labs(tag = "C)") +
  guides(linetype = "none") + theme(legend.position = 'top') +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_size_ctrl)

p_size_R<-ggplot(combo_decarcerate_res %>% filter (mech_num %in% c(4))) +
  geom_line(aes(x=control, y = mean_inf/prob_ob/base_size/PROP_SYMPTOMATIC, col = as.factor(R),linetype = as.factor(mech_num)), size = 1) +
  geom_line(data = simple, aes(x=control, y = 1-control), size = 1, color = 'grey') +
  # geom_line(data = simple, aes(x=control, y = .995-control), linetype = 'dashed', size = 1) +
  geom_point(data = combo_decarcerate_res_sparse %>% filter (mech_num %in% c(4)), aes(x=control, y = mean_inf/prob_ob/base_size/PROP_SYMPTOMATIC, col = as.factor(R), shape = as.factor(R),linetype = as.factor(mech_num)), size = 2.5) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  scale_color_manual(name = 'Baseline R', values = man_color) +
  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized outbreak size') + labs(tag = "D)") +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
#(p_size_R)

g_component <- grid.arrange(p_prob_ctrl,p_size_ctrl,p_prob_R, p_size_R, nrow = 2)
ggsave('Figs/depop_component.jpg',plot = g_component, scale = 1.25)
ggsave("Figs/TIFF/Figure 3.eps",plot = g_component, width=4320, units="px", dpi=300)


###############################
# Same sample numbers
###############################

control_0 <- combo_decarcerate_res %>% filter (control == 0, mech_num == 3)
control_20 <- combo_decarcerate_res %>% filter (control == 0.2, mech_num == 3)
1- control_20$mean_inf/control_0$mean_inf

###############################
# Old stuff
###############################

# p_tot<-ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = mean_inf, col = as.factor(R), linetype = as.factor(mech_num))) +
#   geom_point(data = combo_decarcerate_res_sparse, aes(x=control, y = mean_inf, col = as.factor(R), shape = as.factor(R), linetype = as.factor(mech_num))) +
#   #  geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
#   scale_linetype_discrete(name = NULL, limits = as_factor(MECHANISM_TABLE$mech_num), labels = MECHANISM_TABLE$`Mechanism of control`) +
#   xlab(expression('Degree of control,' ~ gamma)) + ylab('Infections expected\nper 100 days') +
#   guides(color = "none", shape = "none") + theme(legend.position = "top") +
#   theme(legend.text = element_text(size = 14),legend.title = element_text(size = 16),text = element_text(size=20))
# (p_tot)
# 
# 
# p_tot_n<-
#   ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = mean_inf/base_inf, col = as.factor(R),linetype = as.factor(mech_num))) +
#   geom_point(data = combo_decarcerate_res_sparse, aes(x=control, y = mean_inf/base_inf, col = as.factor(R), shape = as.factor(R),linetype = as.factor(mech_num))) +
#   geom_line(data = simple, aes(x=control, y = 1-control)) +
#   scale_color_discrete(name = 'Baseline R') +
#   scale_shape_discrete(name = 'Baseline R') +
#   xlab(expression('Degree of control,' ~ gamma)) + ylab('Normalized\ninfection rate') +
#   guides(linetype = "none") + theme(legend.position = 'top') +
#   theme(legend.text = element_text(size = 10),legend.title = element_text(size = 16),text = element_text(size=20))
# (p_tot_n)
# 
# g_impact_old <- grid.arrange(p_tot,p_tot_n, nrow = 2)
# ggsave('Figs/depop_impact_oldstyle.jpg',plot = g_impact_old, scale = 1.2)

###############################

# p_prob<-ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = prob_ob, col = as.factor(R), linetype = as.factor(prop_factor))) +
#   geom_point(data = combo_decarcerate_res_sparse, aes(x=control, y = prob_ob, col = as.factor(R), shape = as.factor(R), linetype = as.factor(prop_factor))) +
#   #  geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
#   scale_color_discrete(name = 'Baseline R') +
#   scale_shape_discrete(name = 'Baseline R') +
#   #  scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
#   xlab('Reduction of susceptible population') + ylab('Probabilty of an an outbreak') +
#   guides(linetype = "none") + theme(legend.position = 'top') +
#   theme(legend.title = element_text(size = 16),text = element_text(size=20))
# (p_prob)
# 
# p_prob_n<-ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = prob_ob/base_prob, col = as.factor(R),linetype = as.factor(prop_factor))) +
#   geom_line(data = simple, aes(x=control, y = 1-control)) +
#   geom_point(data = combo_decarcerate_res_sparse, aes(x=control, y = prob_ob/base_prob, col = as.factor(R),  shape = as.factor(R),linetype = as.factor(prop_factor))) +
#   #  scale_color_discrete(name = 'Baseline R') +
#   scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
#   xlab('Reduction of susceptible population') + ylab('Normalized outbreak probability') +
#   theme(legend.position = "none") +
#   theme(legend.title = element_text(size = 16),text = element_text(size=20))
# (p_prob_n)
# 
# p_size<-ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = ob_size, col = as.factor(R), linetype = as.factor(prop_factor))) +
#   geom_point(data = combo_decarcerate_res_sparse, aes(x=control, y = ob_size, col = as.factor(R), shape  = as.factor(R), linetype = as.factor(prop_factor))) +
#   #  geom_line(data = simple, aes(x=control, y = mean_inf, col = as.factor(R)), linetype = "dotted") +
#   scale_color_discrete(name = 'Baseline R') +
#   scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
#   xlab('Reduction of susceptible population') + ylab('Outbreak size') +
#   guides(color = "none", shape = "none") + theme(legend.position = 'top') +
#   theme(legend.title = element_text(size = 16),text = element_text(size=20))
# (p_size)
# 
# p_size_n<-ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = ob_size/base_size, col = as.factor(R),linetype = as.factor(prop_factor))) +
#   geom_line(data = simple, aes(x=control, y = 1-control)) +
#   geom_point(data = combo_decarcerate_res_sparse %>% filter(prop_factor == 1) , aes(x=control, y = ob_size/base_size, col = as.factor(R), shape = as.factor(R),linetype = as.factor(prop_factor))) +
#   scale_color_discrete(name = 'Baseline R') +
#   scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
#   xlab('Reduction of susceptible population') + ylab('Normalized outbreak size') +
#   theme(legend.position = "none") +
#   theme(legend.title = element_text(size = 16),text = element_text(size=20))
# (p_size_n)
# 
# g_component <- grid.arrange(p_prob,p_size,p_prob_n,p_size_n, nrow = 2)
# ggsave('Figs/depop_component.jpg',plot = g_component, scale = 1.25)
# 
###############################
# For MIRA
###############################
# 
# p_tot_MIRA<-ggplot(combo_decarcerate_res) +
#   geom_line(aes(x=control, y = mean_inf/base_inf, col = as.factor(R),linetype = as.factor(prop_factor)), size = 2.5) +
#   geom_line(data = simple, aes(x=control, y = 1-control), size = 2.5) +
#   scale_color_discrete(name = 'Baseline R') +
#   guides(linetype = FALSE) +
#   scale_linetype_discrete(name = 'R variable', limits = as_factor(c(1,0)), labels = c('Yes','No')) +
#   xlab('Reduction of susceptible population') + ylab('Normalized infection rate') +
# #  theme(legend.position = "top") +
#   theme(legend.title = element_text(size = 16),text = element_text(size=24))
# (p_tot_MIRA)
# 
# #g_impact <- grid.arrange(p_prob,p_size,p_tot,p_prob_n,p_size_n,p_tot_n, nrow = 2)
# ggsave('Figs/depop_impact_MIRA.jpg',plot = p_tot_MIRA, scale = 0.7)


