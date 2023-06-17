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
### Parameters' possible values:
### n
N_node_persubj_list = list(30)
### beta
conn_patt_sep_list = list(1.9)


top_level_folder = "../Results/Rdata"
setup = 'rerun_Initialization'
default_setting = 'pr=0.9,n=30,beta=1.9,V=80'

method = 'main_v5_cdf'
for (id_split in 1:split) {
  for (N_restart in c(1,2,3,10)) {
    ### beta
    for (i in 1:length(conn_patt_sep_list)) {
      conn_patt_sep = conn_patt_sep_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = id_split*100+j+20000
        tryCatch(main_v5_cdf(SEED = SEED,
                             N_node_vec = rep(30,1),
                             conn_prob_mean = 0.9,
                             conn_patt_sep = 1.9,
                             time_shift_mean_vec = rep(40,3),
                             t_vec = seq(0,200,length.out=200),
                             freq_trun_vec=c(Inf),
                             gamma = 0.01,
                             save_est_history = TRUE,
                             rand_init=TRUE,
                             N_restart = N_restart,
                             conv_thres = 0, MaxIter = 10,
                             N_clus_min = 3, N_clus_max = 3),
                 error = function(x) print(SEED))
      }
      param_name = "beta"
      param_value = conn_patt_sep
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', 'N_restart_v2', '/', N_restart,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
  }
  ### beta
  for (i in 1:length(conn_patt_sep_list)) {
    conn_patt_sep = conn_patt_sep_list[[i]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = id_split*100+j+20000
      tryCatch(main_v5_cdf(SEED = SEED,
                           N_node_vec = rep(30,1),
                           conn_prob_mean = 0.9,
                           conn_patt_sep = 1.9,
                           time_shift_mean_vec = rep(40,3),
                           t_vec = seq(0,200,length.out=200),
                           freq_trun_vec=c(Inf),
                           gamma = 0.01,
                           save_est_history = TRUE,
                           our_restart = TRUE,
                           # N_restart = N_restart,
                           conv_thres = 0, MaxIter = 10,
                           N_clus_min = 3, N_clus_max = 3),
               error = function(x) print(SEED))
    }
    param_name = "beta"
    param_value = conn_patt_sep
    folder_path = paste0(top_level_folder, '/', setup, '/', method,
                         '/', 'Our_init_v7',
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  
}


