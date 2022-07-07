# ASSUMPTIONS
# Discrete, synchronized generations with deterministic rate of infection
# Within each block, transmission events are independent and identical
# Prisoners do not move from one block to another

library(here)
library(tidyverse)
library (gridExtra)

rm(list = ls())

source(here('0 - prison_fxns.R'))

#############################################################

r2p5 <- gen_time_series(R_max = 2.5)
r1p5 <- gen_time_series(prop_capacity = 1, R_max = 1.5)
r2p5_long <- r2p5 %>% pivot_longer(c(names(r2p5)[-1]), names_to = 'Var') %>% mutate (R = 2.5)
r1p5_long <- r1p5 %>% pivot_longer(c(names(r2p5)[-1]), names_to = 'Var') %>% mutate (R = 1.5)

matcher <- tibble(Var = c('num_exp','num_inf','num_susc','num_rec'), compartment = c('Exposed','Infectious','Susceptible','Removed'))
temp_data <- bind_rows(r2p5_long,r1p5_long) %>% filter(Var != 'num_vac') %>% right_join(matcher)

time_pts <- unique(temp_data$time)
sparse_time <- time_pts[seq(20,length(time_pts),40)]
temp_data2 <- temp_data %>% filter(time %in% sparse_time)

ggplot(data = temp_data) +
  geom_line(aes(x=time, y = value, col = compartment)) +
  geom_point(data = temp_data2, aes(x=time, y = value, col = compartment, shape = compartment)) +
  facet_wrap(~R) +
  xlab("Time (days)") + ylab("Population size") +
  labs(col = "Infectious\nstate", shape = "Infectious\nstate") +
  theme(legend.title = element_text(size = 20),text = element_text(size=20)) 
ggsave('Figs/seir_example.jpg')

#decarceration_results <- explore_decarceration(prop_capacity_arr = seq(1,.5, by = -.01), R_max_arr = c(1.5,2,3,5,8))
#write.csv(decarceration_results,"Results/decarceration_results.csv",row.names = FALSE)
decarceration_results <- read_csv("Results/decarceration_results.csv")

# Tot # of hospitalizatin
# Peak of bed use
# Timing of peak
# Peak # of beds

# p1 <- ggplot(data = decarceration_results) +
#   geom_line(aes(x=1-prop_capacity, y = tot_inf, col = as.factor(R_max))) +
#   geom_point(aes(x=1-prop_capacity, y = tot_inf, col = as.factor(R_max), shape = as.factor(R_max))) +
#   ylab ("Total infected") +
#   xlab ("Reduction of susceptible population") +
#   theme(legend.position = "none") + theme(legend.title = element_text(size = 16),text = element_text(size=16)) 
# (p1)

# p2 <- ggplot(data = decarceration_results) +
#   geom_line(aes(x=1-prop_capacity, y = prop_inf, col = as.factor(R_max))) +
#   geom_point(aes(x=1-prop_capacity, y = prop_inf, col = as.factor(R_max), shape = as.factor(R_max))) +
#   ylab ("Proportion infected") +
#   xlab ("Reduction of susceptible population") +
#   theme(legend.position = "none") + theme(legend.title = element_text(size = 16),text = element_text(size=16)) 
# (p2)

p3 <- ggplot(data = decarceration_results) +
  geom_line(aes(x=1-prop_capacity, y = time_inf_peak, col = as.factor(R_max))) +
  geom_point(aes(x=1-prop_capacity, y = time_inf_peak, col = as.factor(R_max), shape = as.factor(R_max))) +
  ylab ("Time of peak\ninfections (days)") +
  xlab(expression('Degree of control,' ~ gamma))  +
  labs(tag = "A)") +
  theme(legend.position = "none") + theme(legend.title = element_text(size = 16),text = element_text(size=16)) 
# (p3)

p4 <- ggplot(data = decarceration_results) +
  geom_line(aes(x=1-prop_capacity, y = max_num_inf, col = as.factor(R_max))) +
  geom_point(aes(x=1-prop_capacity, y = max_num_inf, col = as.factor(R_max), shape = as.factor(R_max))) +
  ylab ("Maximum infected\nat once") +
  xlab(expression('Degree of control,' ~ gamma))  +
  labs(col = "Baseline R", shape = "Baseline R", tag = "B)") + 
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10), text = element_text(size=16)) +
  theme(legend.position = c(0.8, 0.9), legend.direction = 'horizontal')
# (p4)

p_all <- grid.arrange(p3,p4)
ggsave("Figs/seir_res.jpg",plot = p_all)
ggsave("Figs/TIFF/Figure 4.eps",plot = p_all, width=2880, units="px", dpi=300)