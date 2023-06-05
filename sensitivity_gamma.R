#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

# Load libraries ----------------------------------------------------------

library("optparse") 
library(foreach)
library(doParallel)


# User input setup --------------------------------------------------------

option_list = list(
  make_option(c("-n", "--N_trial"), type="integer", default=200, 
              help="number of repeated trials"),
  make_option("--split", type="integer", default=20)
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

N_trial_total = opt$N_trial
split = opt$split


N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
N_clus = 3

### ICL vs gamma, CDF, SAME scale across clusters -----
if (FALSE){
  ### Parameters' possible values:
  ### n
  N_node_persubj_list = list(30)
  ### beta
  conn_patt_sep_list = list(1.9)
  ### gamma
  gamma_value_list = list(0.001,0.003,0.01,0.03,0.1,0.3,1)
  
  top_level_folder = "../Results/Rdata"
  setup = 'sensitivity_gamma'
  method = 'main_v5_cdf_v1'
  
  default_setting = 'pr=0.9,n=30,beta=1.9,V=80'
  for (. in 1:split) {
    ### gamma
    for (i in 1:length(gamma_value_list)) {
      gamma = gamma_value_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_cdf(SEED = SEED,
                             N_node_vec = rep(30,1),
                             conn_prob_mean = 0.9,
                             conn_patt_sep = 1.9,
                             time_shift_mean_vec = rep(40,N_clus),
                             t_vec = seq(0,200,length.out=200),
                             gamma = gamma,
                             freq_trun_vec = c(Inf),
                             MaxIter = 10,
                             N_clus_min = 3,
                             N_clus_max = 3),
                 error = function(x) print(SEED))
      }
      param_name = "gamma"
      param_value = gamma
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
  }
  
}


### ICL vs gamma, CDF, DIFFERENT scale across clusters -----

### Parameters' possible values:
### n
N_node_persubj_list = list(30)
### beta
conn_patt_sep_list = list(1.9)
### gamma
gamma_value_list = list(0.001,0.003,0.01,0.03,0.1,0.3,1)

top_level_folder = "../Results/Rdata"
setup = 'sensitivity_gamma'
method = 'main_v5_cdf_v1'

default_setting = 'pr=vary,n=30,beta=1.9,V=80'
for (. in 1:split) {
  ### gamma
  for (i in 1:length(gamma_value_list)) {
    gamma = gamma_value_list[[i]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_cdf(SEED = SEED,
                           N_node_vec = rep(30,1),
                           conn_prob_mean = 0.5,
                           conn_prob_rad = 0.4,
                           conn_patt_sep = 1.9,
                           time_shift_mean_vec = rep(40,N_clus),
                           t_vec = seq(0,200,length.out=200),
                           gamma = gamma,
                           freq_trun_vec = c(Inf),
                           MaxIter = 10,
                           N_clus_min = 3,
                           N_clus_max = 3),
               error = function(x) print(SEED))
    }
    param_name = "gamma"
    param_value = gamma
    folder_path = paste0(top_level_folder, '/', setup, '/', method,
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
}







