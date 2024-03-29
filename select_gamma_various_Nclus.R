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
  make_option(c("-n", "--N_trial"), type="integer", default=1000, 
              help="number of repeated trials"),
  make_option("--split", type="integer", default=100)
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


### Parameters' possible values:
### gamma
gamma_value_list = list(0.001,0.003,0.01,0.03,0.1,0.3,1,3,10)

top_level_folder = "../Results/Rdata"
setup = 'select_gamma_various_Nclus_v2'

for (id_split in 1:split){
  # ICL vs gamma, CDF, IDENTICAL scale across clusters -----
  if (TRUE) {
    default_setting = 'pr=0.9,n=30,beta=1.9,V=80'
    for (N_clus_est in c(2,3,4)){
      method = paste0('main_v5_cdf_Nclusest', N_clus_est)
      for (i in 1:length(gamma_value_list)) {
        gamma = gamma_value_list[[i]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = id_split*100+j+200000
          print(SEED)
          tryCatch(main_v5_cdf(SEED = SEED,
                               N_node_vec = rep(30,1),
                               conn_prob_mean = 0.9,
                               conn_patt_sep = 1.9,
                               time_shift_mean_vec = rep(40,N_clus),
                               t_vec = seq(0,200,length.out=200),
                               gamma = gamma,
                               freq_trun_vec = c(Inf),
                               MaxIter = 10,
                               N_clus_min = N_clus_est,
                               N_clus_max = N_clus_est),
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
  # ICL vs gamma, CDF, DIFFERENT scale and DIFFERENT shape across clusters -----
  if (TRUE) {
    default_setting = 'pr=vary,n=30,beta=1.9,V=80'
    for (N_clus_est in c(2,3,4)){
      method = paste0('main_v5_cdf_Nclusest', N_clus_est)
      for (i in 1:length(gamma_value_list)) {
        gamma = gamma_value_list[[i]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = id_split*100+j+200000
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
                               N_clus_min = N_clus_est,
                               N_clus_max = N_clus_est),
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
  # ICL vs gamma, CDF, DIFFERENT scale and IDENTICAL shape across clusters -----
  default_setting = 'pr=vary,n=30,beta=1,V=80'
  for (N_clus_est in c(2,3,4)){
    method = paste0('main_v5_cdf_Nclusest', N_clus_est)
    for (i in 1:length(gamma_value_list)) {
      gamma = gamma_value_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = id_split*100+j+200000
        tryCatch(main_v5_cdf(SEED = SEED,
                             N_node_vec = rep(30,1),
                             conn_prob_mean = 0.5,
                             conn_prob_rad = 0.4,
                             conn_patt_sep = 1,
                             time_shift_mean_vec = rep(40,N_clus),
                             t_vec = seq(0,200,length.out=200),
                             gamma = gamma,
                             freq_trun_vec = c(Inf),
                             MaxIter = 10,
                             N_clus_min = N_clus_est,
                             N_clus_max = N_clus_est),
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

