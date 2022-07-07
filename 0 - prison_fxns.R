library(tidyverse)
library (gridExtra)

# Regarding Number introductions / day
PREV <- c(1e-4)
N_C <- 10 # staff contacts / per resident
SAR <- 0.01
N_0 <- 1000

# Probability of an outbeak occuring
MAX_CONTAINED <- 10
INTRO_NUM <- 1
K_S <- .2

# Dynamics of outbreak
NUM_INDEX <- 1
TIME_STEP <- .2
INC_TIME <- 3
INF_TIME <- 7
PROP_SYMPTOMATIC <- 0.75
SIM_TIME <- 100 # Days allowed for an introduction to occur
VAC_TRIGGER <- 10

################

#Simulate clusters that occur after an introduction
#Outputs is a table theat lists the number of clusters of an individual size, and the numeber that lead to full blow outbreaks (last row of output)
#Inputs are:
#   the reproduciton number (R)
#   the dispersion parameter (k)
#   the size at which point an outbreak is deemed to occur (thresh)
#   the number of introductions to simulate (number)
sim_cluster <-function(R,k,thresh,number) {
  all_sizes <- tibble(i = 1:number) %>% group_by(i) %>% do({
    num_inf <- 1
    cluster_size <- 1
    while (num_inf > 0) {
      num_inf <- sum(rnbinom(num_inf,size = k, mu = R))
      cluster_size <- cluster_size + num_inf
      if(cluster_size > thresh) {
        cluster_size <- thresh + 1e10
        num_inf <- 0
      }
    }
    tibble(cluster_size = cluster_size)
  })
  all_sizes %>% group_by(cluster_size) %>% summarize(n = n())
}

# prob that i cases cause j cases
calc_r_ij <- function (i,j,r,k) {
  log_r_ij <- lgamma(j +k*i) - lgamma(j+1) - lgamma(k*i) + 
    k * i * log(k/(r+k)) +
    j * log(r/(r+k))
  r_ij<- exp(log_r_ij)
}

# prob a single introduction results in an extinction
calc_prob_ext <- function(thresh, R, k) {
  clust_size=1:thresh
  s_ij_arr <- calc_r_ij(i=clust_size,j=clust_size-1,R,k)/clust_size
  p_extinct <- sum(s_ij_arr)
}

# Likelihood of data
#
# thresh denotes the size cutoff that determines when a cluster becomes an outbreak.
# c_j_arr = probability of having a cluster of size j
calc_cluster_logL <- function(data, thresh, R, k) {
  cluster_size <- 1:thresh
  c_j_arr <- calc_r_ij(i=cluster_size,j=cluster_size-1,R,k)/cluster_size
  cluster_size[thresh+1] <- thresh + 1
  c_j_arr[thresh+1] <- 1-sum(c_j_arr)
  prob_arr <- tibble(cluster_size = cluster_size, prob = c_j_arr)
  
  # Truncate right tail
  data$cluster_size <- ifelse(data$cluster_size > thresh,thresh+1, data$cluster_size)

  #Do log L
  data <- left_join(data,prob_arr,by = 'cluster_size')
  logL <- sum(data$n * log(data$prob))
}

gen_time_series <- function(prop_capacity = 1, unit_size = N_0, R_max, prop_factor = 1, trigger = VAC_TRIGGER, vac_delay = 0, ve = 0, time_step = TIME_STEP) {
  
  current_time <-  tibble(time = 0, num_exp = NUM_INDEX,  num_inf = 0, num_susc = unit_size*prop_capacity - NUM_INDEX, num_rec = 0, num_vac = 0) #stores results for current generation
  time_series <- current_time #Will store all data
  
  R <- R_max * (1 - (1 - prop_capacity) * prop_factor)
  N <- unit_size * prop_capacity
  
  vac_status <- 0
  num_vac <- 0
  vac_time <- Inf
  while (current_time$num_exp > NUM_INDEX/1e3) {
    last_time <- current_time
    if(trigger > 0 && unit_size*prop_capacity - last_time$num_susc >= trigger && is.infinite(vac_time)) {
      vac_time <-last_time$time + vac_delay
    }
    if(last_time$time + TIME_STEP > vac_time && vac_status == 0) {
      #    if(unit_size*prop_capacity - last_time$num_susc > trigger && vac_status == 0) {
      vac_status <- 1
      num_vac <- last_time$num_susc * ve
      last_time$num_susc <- last_time$num_susc - num_vac
    }
    new_exp <- R*(time_step)*(last_time$num_inf/INF_TIME) * last_time$num_susc / N
    new_inf <- last_time$num_exp * (time_step/INC_TIME)

    num_susc <- last_time$num_susc - new_exp # Susceptibles are decreased by the # who are now infectious or vaccinated
    num_exp <- last_time$num_exp + new_exp - new_inf
    num_inf <- last_time$num_inf + new_inf - last_time$num_inf * (time_step/INF_TIME)
    num_rec <- last_time$num_rec + last_time$num_inf * (time_step/INF_TIME)
    
    current_time <- tibble(time = last_time$time + time_step, num_susc = num_susc, num_exp = num_exp,  
                           num_inf = num_inf, num_rec = num_rec, num_vac = num_vac)
    time_series <- bind_rows(time_series,current_time)
  }
  time_series
}

explore_decarceration <- function(prop_capacity_arr, R_max_arr, prop_factor = 1, unit_size = N_0, time_step = TIME_STEP) {
  setup_arr <-tibble(prop_capacity = rep(prop_capacity_arr, time = length(R_max_arr)), R_max = rep(R_max_arr, each = length(prop_capacity_arr))) 
  results <- setup_arr %>% group_by(prop_capacity,R_max) %>% do({
    time_series <- gen_time_series(prop_capacity = .$prop_capacity, R_max = .$R_max, prop_factor = prop_factor, unit_size = unit_size, time_step = time_step)
    tot_inf <- last(time_series$num_rec)
    prop_inf <- tot_inf / unit_size / .$prop_capacity
    time_inf_peak <- time_series$time[which.max(time_series$num_inf)]
    max_num_inf <- max(time_series$num_inf)
    
    tibble(tot_inf = tot_inf, prop_inf = prop_inf, time_inf_peak = time_inf_peak, max_num_inf = max_num_inf)
  })
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}