library(here)
rm(list=ls())
source(here('0 - prison_fxns.R'))

R0_ARR <- c(1.5, 2.5, 3.5)
VAC_EFF_ARR <- c(.2, .5, .8)
VAC_DELAY_ARR <- seq(5,50,5)
TRIGGER_ARR <- c(2,5,10)

eval_vaccination <- function (vac_delay_arr, R0_arr, vac_eff_arr, trigger_arr) {
  tibble (expand_grid(vac_delay = vac_delay_arr, R0 = R0_arr, vac_eff = vac_eff_arr, trigger = trigger_arr)) %>% group_by(vac_delay, R0, vac_eff, trigger) %>% do ({

    num_ic <- SIM_TIME * PREV * N_0 * N_C * SAR #Number of infectious contacts
    R0 <- .$R0[1]
    prob_ext <- calc_prob_ext(MAX_CONTAINED, R0, K_S)
    #Combined
    prob_ob <- 1 - prob_ext^num_ic
    
    time_series <- gen_time_series(R_max = R0, vac_delay = .$vac_delay[1], ve = .$vac_eff[1],trigger = .$trigger[1])
    size_ob <- last(time_series$num_rec)
    tibble(prob_ob = prob_ob, prob_ext = prob_ext, size_ob = size_ob, total_inf = prob_ob*size_ob) 
  })
}

eval_preemp_vaccination <- function (R0_arr, vac_eff_arr) {
  tibble (expand_grid(R0 = R0_arr, vac_eff = vac_eff_arr)) %>% group_by(R0, vac_eff) %>% do ({
    
    num_ic <- SIM_TIME * PREV * N_0 * N_C * SAR *(1-.$vac_eff[1]) #Number of infectious contacts
    R0 <- .$R0[1]
    prob_ext <- calc_prob_ext(MAX_CONTAINED, R0*(1-.$vac_eff[1]), K_S)
    #Combined
    prob_ob <- 1 - prob_ext^num_ic
    
    time_series <- gen_time_series(prop_capacity = (1-.$vac_eff[1]),R_max = R0, vac_delay = 0,trigger = 0)
    size_ob <- last(time_series$num_rec)
    tibble(prob_ob = prob_ob, prob_ext = prob_ext, size_ob = size_ob, total_cases = prob_ob*size_ob*PROP_SYMPTOMATIC) 
  })
}

vac_res <- eval_vaccination(vac_delay_arr = VAC_DELAY_ARR,R0_arr = R0_ARR, vac_eff_arr = VAC_EFF_ARR, trigger_arr =TRIGGER_ARR)
write.csv(vac_res,"Results/vac_res.csv",row.names = FALSE)
vac_res <- read_csv("Results/vac_res.csv")
vac_res <- vac_res %>% mutate(total_cases = total_inf*PROP_SYMPTOMATIC)

# ggplot(vac_res %>% filter(trigger == 10)) +
#   geom_line(aes(x=vac_delay,y=total_inf,col=as_factor(R0), linetype = as_factor(vac_eff)))

PREEMP_VAC_ARR <- seq(0,0.2,.1)
preemp_vac_res <- eval_preemp_vaccination(R0_arr = R0_ARR, vac_eff_arr = PREEMP_VAC_ARR)
preemp_vac_res <- bind_rows(preemp_vac_res %>% mutate(vac_delay = min(VAC_DELAY_ARR)),preemp_vac_res %>% mutate(vac_delay = max(VAC_DELAY_ARR)))

ggplot(vac_res %>% filter(trigger == 10)) +
  geom_line(aes(x=vac_delay,y=total_cases,col=as_factor(vac_eff)), size = 1) +
  geom_point(aes(x=vac_delay,y=total_cases,col=as_factor(vac_eff),shape=as_factor(vac_eff)), size = 2.5) +
  geom_line(data = preemp_vac_res, aes(x = vac_delay,y = total_cases,linetype = as_factor(vac_eff)), size = 1) +
  xlab('Delay of implementing reactive control (days)') + ylab('Cases expected per 100 days') +
  scale_linetype_discrete(name = 'Preemptive\nvaccine\nefficacy') +
  scale_color_discrete(name = 'Reactive\nvaccine\nefficacy') +
  scale_shape_discrete(name = 'Reactive\nvaccine\nefficacy') +
  theme(legend.title = element_text(size = 16),text = element_text(size=20)) +
  facet_wrap(~ R0)
ggsave('Figs/reactive_control.jpg', scale = 1.25)
ggsave("Figs/TIFF/Figure 6.eps",width=4320, units="px", dpi=300)