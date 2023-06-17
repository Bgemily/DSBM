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

### ARI vs SNR, V!=0 -----

### Parameters' possible values:
### n
N_node_persubj_list = list(30,42,54,66,78,90)
### beta
conn_patt_sep_list = list(1.3,1.4,1.5,1.6,1.7,1.8,1.9)
### V
time_shift_mean_vec_list = list(rep(15,N_clus), rep(20,N_clus),
                                rep(25,N_clus), rep(30,N_clus),
                                rep(35,N_clus), rep(40,N_clus))


top_level_folder = "../Results/Rdata"
setup = 'rerun_sim_SNR_Vnot0'
default_setting = 'pr=0.9,n=30,beta=1.3,V=80'

for (id_split in 1:split) {
  method = 'main_v5_pdf_Nclusest1_kernel_v2'
  for (freq_trun in c(4)){
    ### N_node
    for (i in 1:length(N_node_persubj_list)) {
      N_node = N_node_persubj_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = id_split*100+j+20000
        tryCatch(main_v5_pdf_Nclusest1(SEED = SEED,
                             N_node_vec = rep(N_node,1),
                             conn_prob_mean = 0.9,
                             conn_patt_sep = 1.3,
                             time_shift_mean_vec = rep(40,N_clus),
                             t_vec = seq(0,200,length.out=200),
                             freq_trun_vec = c(freq_trun),
                             MaxIter = 10,
                             gamma = 0.0001,
                             N_clus_min = 1, N_clus_max = 1),
                 error = function(x) print(SEED))
      }
      param_name = "n"
      param_value = N_node
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           # '/','freq_trun','/',freq_trun,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### beta
    for (i in 1:length(conn_patt_sep_list)) {
      conn_patt_sep = conn_patt_sep_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = id_split*100+j+20000
        tryCatch(main_v5_pdf_Nclusest1(SEED = SEED,
                             N_node_vec = rep(30,1),
                             conn_prob_mean = 0.9,
                             conn_patt_sep = conn_patt_sep,
                             time_shift_mean_vec = rep(40,N_clus),
                             t_vec = seq(0,200,length.out=200),
                             freq_trun_vec = c(freq_trun),
                             MaxIter = 10,
                             gamma = 0.0001,
                             N_clus_min = 1, N_clus_max = 1),
                 error = function(x) print(SEED))
      }
      param_name = "beta"
      param_value = conn_patt_sep
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           # '/','freq_trun','/',freq_trun,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
  }
  
  
}



