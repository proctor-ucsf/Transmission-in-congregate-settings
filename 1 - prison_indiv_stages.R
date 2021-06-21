library(tidyverse)
library (gridExtra)
rm(list=ls())

#############################################################

source('0 - prison_fxns.R')

#############################################################
# Number of introduction per day
#############################################################

N_C_ARR <- seq(2,10,.1) # Staff contacts per resident
SAR_ARR <- seq(.01, .05, .001)
#K_P_ARR <- 1e3#c(0.2,1e3)

res_intro_arr <- tibble(N_c = N_C_ARR) %>% group_by (N_c) %>% do ({
  N_ic <- PREV * N_0 * .$N_c[1]
  tibble(SAR = SAR_ARR) %>% group_by (SAR) %>% do ({
    SAR <- .$SAR[1]
    approx_days <- 1/(N_ic*SAR)
    tibble(approx_days = approx_days)
  })
})

g1a <- ggplot(res_intro_arr) +
  scale_fill_distiller(palette = "Spectral") +
  geom_tile(aes(x = N_c, y = SAR, fill = -log(approx_days))) +
#  scale_fill_gradient(low="red", high="blue") +
  xlab("Contacts per prisoner") + ylab("Attack rate") + theme(text = element_text(size=20)) +
  labs(fill = "Average\ndays until\nintroduction\n") + theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))
(g1a)  

g1b <- ggplot(res_intro_arr %>% filter(SAR %in% seq(0.01, 0.05, 0.01))) +
  geom_point(aes(x = N_c, y = approx_days, col = as.factor(SAR))) +
  xlab("Contacts per prisoner") + ylab("Average days\nuntil introduction") + theme(text = element_text(size=20)) +
  labs(col = "Transmission\nprobability") + theme(legend.title = element_text(size = 12),legend.text = element_text(size = 12))
(g1b)

#g_intro <- grid.arrange(g1a,g1b)
ggsave('Figs/prob_intro.jpg',plot = g1b)

#############################################################
# Probability of an introduction leading to an outbreak
#############################################################
INTRO_ARR <- c(1,2,5)

R_ARR <- seq (1.5, 3, .1)
K_ARR <- c(0.2, .5, 1, 1e3)

res_prob_ob_arr <- tibble(R = R_ARR) %>% group_by (R) %>% do ({
  R <- .$R[1]
  tibble(k = K_ARR) %>% group_by (k) %>% do ({
    p_extinct <- calc_prob_ext(thresh = MAX_CONTAINED, R = R, k =.$k[1])
    tibble(num_intro = INTRO_ARR,p_ob = 1-p_extinct^INTRO_ARR)
  })
})

g_prob_ob <- ggplot(res_prob_ob_arr) +
  geom_line(aes(x=R, y = p_ob, col = as.factor(k))) +
  facet_wrap(~num_intro) +
  xlab("Reproduction number (R)") + ylab("Outbreak\nprobability") +
  #theme(text = element_text(size=20)) +
  labs(col = "Dispersion\nparameter (k)")
(g_prob_ob)
ggsave('Figs/prob_ob.jpg',plot = g_prob_ob)


#############################################################
res_prob_ob_arr <- tibble(R = R_ARR) %>% group_by (R) %>% do ({
  R <- .$R[1]
  tibble(k = K_ARR) %>% group_by (k) %>% do ({
    p_extinct <- calc_prob_ext(thresh = MAX_CONTAINED, R = R, k =.$k[1])
    tibble(p_ob = 1-p_extinct)
  })
})

rename_arr <- tribble(~k, ~kk,
                      0.2, 'Superspreader (0.2)',
                      0.5, 'Heterogeneous (0.5)',
                      1, 'Geometric (1)',
                      1000, 'Homogeneous (Infinite)')

r2 <- right_join(res_prob_ob_arr,rename_arr, by = 'k')

g_prob_ob <- ggplot(r2) +
  geom_line(aes(x=R, y = p_ob, col = as.factor(kk))) +
  xlab("Reproduction number (R)") + ylab("Outbreak\nprobability") +
  #theme(text = element_text(size=20)) +
  labs(col = "Dispersion\nparameter (k)") + theme(legend.title = element_text(size = 18),text = element_text(size=30)) 
(g_prob_ob)

ggsave('Figs/prob_ob_single.jpg',plot = g_prob_ob)

#############################################################
# Size of outbreak
#############################################################
# Size of an outbreak (R, N_p)
source('1A - prison_seir.R')

