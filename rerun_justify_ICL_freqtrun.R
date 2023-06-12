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

### Parameters' possible values:
### beta
conn_patt_sep_list = list(1.9)

top_level_folder = "../Results/Rdata"
setup = 'rerun_justify_ICL_freqtrun'
default_setting = 'pr=0.9,n=90,beta=1.9,V=80'


for (id_split in 1:split) {
  # CDF
  method = 'main_v5_cdf'
  for (i in 1:length(conn_patt_sep_list)) {
    conn_patt_sep = conn_patt_sep_list[[i]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = id_split*10000+j
      tryCatch(main_v5_cdf(SEED = SEED,
                           N_node_vec = rep(90,1),
                           conn_prob_mean = 0.9,
                           conn_patt_sep = conn_patt_sep,
                           time_shift_mean_vec = rep(40,N_clus),
                           t_vec = seq(0,200,length.out=200),
                           freq_trun_vec = c(Inf),
                           MaxIter = 10,
                           gamma = 0.01,
                           N_clus_min = 1,
                           N_clus_max = 5),
               error = function(x) print(SEED))
    }
    param_name = "beta"
    param_value = conn_patt_sep
    folder_path = paste0(top_level_folder, '/', setup, '/', method,
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  # PDF
  for (freq_trun in 2:9) {
    method = paste0('main_v5_pdf_freqtrun',freq_trun)
    for (i in 1:length(conn_patt_sep_list)) {
      conn_patt_sep = conn_patt_sep_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = id_split*10000+j
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node_vec = rep(90,1),
                             conn_prob_mean = 0.9,
                             conn_patt_sep = conn_patt_sep,
                             time_shift_mean_vec = rep(40,N_clus),
                             t_vec = seq(0,200,length.out=200),
                             freq_trun = freq_trun,
                             MaxIter = 10,
                             gamma = 0.0001,
                             N_clus_min = 1, N_clus_max = 5),
                 error = function(x) print(SEED))
      }
      param_name = "beta"
      param_value = conn_patt_sep
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








