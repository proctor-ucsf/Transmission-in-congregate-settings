# ASSUMPTIONS
# Discrete, synchronized generations with deterministic rate of infection
# Within each block, transmission events are independent and identical
# Prisoners do not move from one block to another

library(tidyverse)
library (gridExtra)

rm(list = ls())

source('0 - prison_fxns.R')

#############################################################

r2p5 <- gen_time_series(R_max = 2.5)
r1p5 <- gen_time_series(prop_capacity = 1, R_max = 1.5)
r2p5_long <- r2p5 %>% pivot_longer(c(names(r2p5)[-1]), names_to = 'Var') %>% mutate (R = 2.5)
r1p5_long <- r1p5 %>% pivot_longer(c(names(r2p5)[-1]), names_to = 'Var') %>% mutate (R = 1.5)

matcher <- tibble(Var = c('num_exp','num_inf','num_susc','num_rec'), compartment = c('Exposed','Infectious','Susceptible','Removed'))
temp_data <- bind_rows(r2p5_long,r1p5_long) %>% filter(Var != 'num_vac') %>% right_join(matcher)

ggplot(data = temp_data) +
  geom_line(aes(x=time, y = value, col = compartment)) +
  facet_wrap(~R) +
  xlab("Time (days)") + ylab("Population size") +
  theme(legend.title = element_text(size = 20),text = element_text(size=20)) 
ggsave('Figs/seir_example.jpg')

#decarceration_results <- explore_decarceration(prop_capacity_arr = seq(1,.5, by = -.01), R_max_arr = c(1.5,2,3,5,8))
#write.csv(decarceration_results,"Results/decarceration_results.csv",row.names = FALSE)
decarceration_results <- read_csv("Results/decarceration_results.csv")

# Tot # of hospitalizatin
# Peak of bed use
# Timing of peak
# Peak # of beds

p1 <- ggplot(data = decarceration_results) +
  geom_line(aes(x=1-prop_capacity, y = tot_inf, col = as.factor(R_max))) +
  ylab ("Total infected") +
  xlab ("Control") +
  theme(legend.position = "none") + theme(legend.title = element_text(size = 18),text = element_text(size=18)) 
(p1)

p2 <- ggplot(data = decarceration_results) +
  geom_line(aes(x=1-prop_capacity, y = prop_inf, col = as.factor(R_max))) +
  ylab ("Proportion infected") +
  xlab ("Control") +
  theme(legend.position = "none") + theme(legend.title = element_text(size = 18),text = element_text(size=18)) 
(p2)

p3 <- ggplot(data = decarceration_results) +
  geom_line(aes(x=1-prop_capacity, y = time_inf_peak, col = as.factor(R_max))) +
  ylab ("Time of peak infections (days)") +
  xlab ("Control") +
  theme(legend.position = "none") + theme(legend.title = element_text(size = 18),text = element_text(size=18)) 
(p3)

p4 <- ggplot(data = decarceration_results) +
  geom_line(aes(x=1-prop_capacity, y = max_num_inf, col = as.factor(R_max))) +
  ylab ("Maximum infected at once") +
  xlab ("Control") +
  labs(col = "Reproduction\nnumber") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10), text = element_text(size=18)) +
  theme(legend.position = c(0.80, 0.70))
(p4)

p_all <- grid.arrange(p1,p2,p3,p4)
ggsave("Figs/seir_res.jpg",plot = p_all)
