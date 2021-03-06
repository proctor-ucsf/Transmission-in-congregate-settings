library(here)
rm(list=ls())
source(here('0 - prison_fxns.R'))

eval_proportion <- function (capacity_arr, prop_N1_arr,R) {
  tibble (prop_N1 = prop_N1_arr) %>% group_by(prop_N1) %>% do ({
    CAP_1 <- capacity_arr[1]
    CAP_2 <- capacity_arr[2]
    N1_size <- (CAP_1 + CAP_2)*.$prop_N1[1]
    N2_size <- (CAP_1 + CAP_2)*(1-.$prop_N1[1])
    
    num_ic <- SIM_TIME * PREV * N1_size * N_C * SAR #Number of infectious contacts
    effective_R1 <- R * N1_size/CAP_1
    prob_ext <- calc_prob_ext(MAX_CONTAINED, effective_R1, K_S)
    #Combined
    prob_ob1 <- 1 - prob_ext^num_ic
    
    num_ic <- SIM_TIME * PREV * N2_size * N_C * SAR #Number of infectious contacts
    effective_R2 <- R * N2_size/CAP_2
    prob_ext <- calc_prob_ext(MAX_CONTAINED, effective_R2, K_S)
    #Combined
    prob_ob2 <- 1 - prob_ext^num_ic
    
    ########
    N1_res <- explore_decarceration(prop_capacity_arr = 1, R_max_arr = effective_R1, prop_factor = 1, unit_size =N1_size, time_step = TIME_STEP)
    N2_res <- explore_decarceration(prop_capacity_arr = 1, R_max_arr = effective_R2, prop_factor = 1, unit_size =N2_size, time_step = TIME_STEP)
    total_cases <- PROP_SYMPTOMATIC*(prob_ob1*N1_res$tot_inf + prob_ob2*N2_res$tot_inf)
    
    tibble(N1 = N1_size, N2 = N2_size, 
           prob_1 = prob_ob1, prob_2 = prob_ob2, 
           size_1 = N1_res$tot_inf, size_2 = N2_res$tot_inf, 
           cases_1 = PROP_SYMPTOMATIC*prob_ob1*N1_res$tot_inf, cases_2 = PROP_SYMPTOMATIC*prob_ob2*N2_res$tot_inf, 
           prob_joint = prob_1 + prob_2 - prob_1*prob_2, total_cases= total_cases)
  })
}

display_proportion_res <- function(proportion_res, fileroot, sparse_arr) {
  
  res_sparse <- proportion_res %>% filter(prop_N1 %in% sparse_arr)
  pp1 <- ggplot(proportion_res) +
    geom_line(aes(x=prop_N1, y = prob_1, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = prob_1, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = "none") +
    theme(text = element_text(size=20)) +
    xlab('Proportion of residents in Unit A') + ylab('Unit A: Outbreak probability') 
  ps1 <- ggplot(proportion_res) + 
    geom_line(aes(x=prop_N1, y = size_1, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = size_1, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = "none") +
    theme(text = element_text(size=20)) +
    xlab('Proportion of residents in Unit A') + ylab('Unit A: Outbreak size') 
  pt1 <- ggplot(proportion_res) + 
    geom_line(aes(x=prop_N1, y = cases_1, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = cases_1, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = c(0.34,0.79), legend.box = "vertical") +
    theme(text = element_text(size=20), legend.title=element_text(size=16), legend.text=element_text(size=12), legend.spacing.y = unit(0, 'mm')) +
    xlab('Proportion of residents in Unit A') + ylab('Unit A: Cases per 100 days')
  pp2 <- ggplot(proportion_res) + 
    geom_line(aes(x=prop_N1, y = prob_2, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = prob_2, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = "none") +
    theme(text = element_text(size=20)) +
    xlab('Proportion of residents in Unit A') + ylab('Unit B: Outbreak probability')
  ps2 <- ggplot(proportion_res) + 
    geom_line(aes(x=prop_N1, y = size_2, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = size_2, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = "none") +
    theme(text = element_text(size=20)) +
    xlab('Proportion of residents in Unit A') + ylab('Unit B: Outbreak size')
  pt2 <- ggplot(proportion_res) + 
    geom_line(aes(x=prop_N1, y = cases_2, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = cases_2, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = "none") +
    theme(text = element_text(size=20)) +
    xlab('Proportion of residents in Unit A') + ylab('Unit B: Cases per 100 days')
  psum <- grid.arrange(pp1,ps1,pt1,pp2,ps2,pt2, nrow = 2)
  (psum)
  ggsave(paste('Figs/',fileroot,'_component.jpg',sep=""),plot = psum, scale = 1.5)
  
  p_jprob <- ggplot(proportion_res) + 
    geom_line(aes(x=prop_N1, y = prob_joint, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = prob_joint, col = Scenario, shape = Scenario), size = 2.5) +
    xlab('Proportion of residents in Unit A') + ylab('Probability of at least\none outbreak occuring') +
    theme(text = element_text(size=20)) +
    geom_point(data = proportion_res %>% group_by(Scenario) %>% filter(prob_joint == min(prob_joint)),
               aes(x = prop_N1, y = prob_joint))
  p_tot <- ggplot(proportion_res) +
    geom_line(aes(x=prop_N1, y = total_cases, col = Scenario), size = 1) +
    geom_point(data = res_sparse, aes(x=prop_N1, y = total_cases, col = Scenario, shape = Scenario), size = 2.5) +
    theme(legend.position = "none") +
    theme(text = element_text(size=20)) +
    xlab('Proportion of residents in Unit A') + ylab('Expected number of\ncases per 100 days') +
    geom_point(data = proportion_res %>% group_by(Scenario) %>% filter(total_cases == min(total_cases)),
               aes(x = prop_N1, y = total_cases))
  psum <- grid.arrange(p_jprob,p_tot)
  ggsave(paste('Figs/',fileroot,'_summary.jpg',sep=""),plot = psum, scale = 1.25)
  (psum)
  # ggsave("Figs/TIFF/Figure 5.eps",plot = psum, width=4320, units="px", dpi=300)
}

scen_1 <- eval_proportion(capacity_arr = c(1000,1000),prop_N1_arr = seq (0.1, 0.9, 0.01), R = 3)
scen_2 <- eval_proportion(capacity_arr = c(500,1500),prop_N1_arr = seq (0.1, 0.9, 0.01), R = 3)
scen_3 <- eval_proportion(capacity_arr = c(800,800),prop_N1_arr = seq (0.1, 0.9, 0.01), R = 2.4)
scen_4 <- eval_proportion(capacity_arr = c(400,400),prop_N1_arr = seq (0.1, 0.9, 0.01), R = 0.75)
scen_5 <- eval_proportion(capacity_arr = c(400,1200),prop_N1_arr = seq (0.1, 0.9, 0.01), R = 2.4)
temp <- SAR
SAR <- 0.05
scen_6 <- eval_proportion(capacity_arr = c(1000,1000),prop_N1_arr = seq (0.1, 0.9, 0.01), R = 6)
SAR <- temp

scen_all <- bind_rows(
  scen_1 %>% mutate(Scenario = "1: R1 = 3, R2 = 3"),
  scen_2 %>% mutate(Scenario = "2: R1 = 6, R2 = 2"),
  scen_3 %>% mutate(Scenario = "3: Scenario 1, with 20% immunity"),
  scen_4 %>% mutate(Scenario = "4: Scenario 1, with 75% immunity"),
  scen_5 %>% mutate(Scenario = "5: Scenario 2, with 20% immunity"),
  scen_6 %>% mutate(Scenario = "6: R1 = 6, R2 = 6, SAR = 0.05")
)
write.csv(scen_all,"Results/scen_all.csv",row.names = FALSE)
scen_all <- read_csv("Results/scen_all.csv")

display_proportion_res(scen_all,'scenario_all',sparse_arr = seq(.15,.85,.1))
