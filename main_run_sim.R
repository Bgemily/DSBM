#!/usr/bin/env Rscript

### Based on v2
### Save results with structure: setup/method/default_setting/param_names/param_values/time_stamp.Rdata
### For example, setup == "SNR_Vis0", "SNR_Vnot0", "Unbalanced_clusters", "Multi_subj"
### method == "our_v3.4"
### default_setting == 'p=0.5,n=30,beta=1.3'
### param_names == 'N_node','beta', 'V'
### param_values == '30','1.3','90'

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


# Simulation setup --------------------------------------------------------

### Common parameters
N_clus = 3
time_shift_struc = max
total_time = 200
max_iter = 3 ### number of iterations when updating time shifts


### Simulation-specific parameters
N_subj = 3 
N_subj_list = list(1,2,3,4,5)

jitter_time_rad_list = list(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)


### n
N_node_persubj_list = list(30,45,60,90,120)
### beta
conn_patt_sep_list = list(1.4,1.5,1.6,1.7,1.8)
### unbalancedness
n_star_list = list(22,18,14,10)
### V
time_shift_mean_vec_list = list(rep(15,N_clus), rep(20,N_clus), 
                                rep(25,N_clus), rep(30,N_clus), 
                                rep(35,N_clus), rep(40,N_clus))

### alpha
conn_patt_var_list = list(1,2,3,4,5,6,7,8)
### p
conn_prob_mean_list = list(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)


# Run simulations ---------------------------------------------------------

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
setup = 'SNR_Vnot0_v4'
method = 'main_v5_pdf_v1'
default_setting = 'pr=0.9,n=30,beta=1.3,V=80'

for (freq_trun in c(7)){
  for (. in 1:split) {
    ### N_node
    for (i in 1:length(N_node_persubj_list)) {
      N_node = N_node_persubj_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5(SEED = SEED,
                         N_node_vec = rep(N_node,1),
                         conn_prob_mean = 0.9,
                         conn_patt_sep = 1.3,
                         time_shift_mean_vec = rep(40,N_clus),
                         t_vec = seq(0,200,length.out=200),
                         freq_trun_vec = c(freq_trun), 
                         N_clus_min = N_clus, N_clus_max = N_clus),
                 error = function(x) print(SEED))
      }
      param_name = "n"
      param_value = N_node
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/','freq_trun','/',freq_trun,
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
        SEED = sample(1:1e7,1)
        tryCatch(main_v5(SEED = SEED,
                         N_node_vec = rep(30,1),
                         conn_prob_mean = 0.9,
                         conn_patt_sep = conn_patt_sep,
                         time_shift_mean_vec = rep(40,N_clus),
                         t_vec = seq(0,200,length.out=200),
                         freq_trun_vec = c(freq_trun), 
                         N_clus_min = N_clus, N_clus_max = N_clus),
                 error = function(x) print(SEED))
      }
      param_name = "beta"
      param_value = conn_patt_sep
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/','freq_trun','/',freq_trun,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### V: time_shift_mean_vec_list
    for (i in 1:length(time_shift_mean_vec_list)) {
      time_shift_mean_vec = time_shift_mean_vec_list[[i]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5(SEED = SEED,
                         N_node_vec = rep(30,1),
                         conn_prob_mean = 0.9,
                         conn_patt_sep = 1.3,
                         time_shift_mean_vec = time_shift_mean_vec,
                         t_vec = seq(0,200,length.out=200),
                         freq_trun_vec = c(freq_trun), 
                         N_clus_min = N_clus, N_clus_max = N_clus),
                 error = function(x) print(SEED))
      }
      param_name = "V"
      param_value = time_shift_mean_vec[1]*2
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/','freq_trun','/',freq_trun,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    
  }
}






# ### ICL vs N_clus and N_basis ----- 
# method = 'main_v5_pdf_v9_multi_Nclus_Nbasis'
# default_setting = 'pr=1,n=90,beta=1.2'
# N_node_persubj_list = list(90)
# 
# for (. in 1:split) {
#   ### N_node
#   for (i in 1:length(N_node_persubj_list)) {
#     N_node = N_node_persubj_list[[i]]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v5(SEED = SEED,
#                        N_node_vec = rep(N_node,1),
#                        conn_prob_mean = 1,
#                        conn_patt_sep = 1.2,
#                        time_shift_mean_vec = rep(20,N_clus),
#                        t_vec = seq(0,200,length.out=200),
#                        freq_trun_vec = c(3,5,7,9), step_size=0.5,
#                        N_clus_min = 1, N_clus_max = 5),
#                error = function(x) print(SEED))
#     }
#     param_name = "n"
#     param_value = N_node
#     folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                          '/', default_setting,
#                          '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#   }
# }
# 
# default_setting = 'pr=1,n=54,beta=1.2'
# N_node_persubj_list = list(54)
# 
# for (. in 1:split) {
#   ### N_node
#   for (i in 1:length(N_node_persubj_list)) {
#     N_node = N_node_persubj_list[[i]]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v5(SEED = SEED,
#                        N_node_vec = rep(N_node,1),
#                        conn_prob_mean = 1,
#                        conn_patt_sep = 1.2,
#                        time_shift_mean_vec = rep(20,N_clus),
#                        t_vec = seq(0,200,length.out=200),
#                        freq_trun_vec = c(3,5,7,9), step_size=0.5,
#                        N_clus_min = 1, N_clus_max = 5),
#                error = function(x) print(SEED))
#     }
#     param_name = "n"
#     param_value = N_node
#     folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                          '/', default_setting,
#                          '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#   }
# }
# 
# 
# ### ARI vs Jitter_radius, V==0 -----
# 
# ### Parameters' possible values: 
# ### n
# # N_node_persubj_list = list(30,42,54,66,78,90)
# N_node_persubj_list = list(90)
# 
# top_level_folder = "../Results/Rdata"
# setup = 'SNR_Vis0'
# method = 'main_v5_pdf_timeshift_jitter_v2'
# default_setting = 'pr=1,n=30,beta=1.15'
# 
# for (jitter_time_rad in c(0,5,10,15,20)) {
#   
#   for (. in 1:split) {
#     ### N_node
#     for (i in 1:length(N_node_persubj_list)) {
#       N_node = N_node_persubj_list[[i]]
#       N_node = 90
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5(SEED = SEED, 
#                          N_node_vec = rep(N_node,1),
#                          conn_prob_mean = 1, 
#                          conn_patt_sep = 1.15,
#                          time_shift_mean_vec = rep(0,N_clus),
#                          t_vec = seq(0,200,length.out=200),
#                          freq_trun=7,
#                          jitter_time_rad = jitter_time_rad,
#                          N_clus_min = N_clus, N_clus_max = N_clus),
#                  error = function(x) print(SEED))
#       }
#       param_name = "jitter_time_rad"
#       param_value = jitter_time_rad
#       folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                            '/', default_setting, 
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     
#   }
#   
# }
# 
# 
# 
# 
# ### N_clus_est, V==0 -----
# 
# ### Parameters' possible values: 
# ### n
# N_node_persubj_list = list(30,42,54,66,78,90)
# # N_node_persubj_list = list(90)
# ### beta
# conn_patt_sep_list = list(1.3,1.4,1.5,1.6,1.7,1.8)
# conn_patt_sep_list = list(1.8)
# 
# top_level_folder = "../Results/Rdata"
# setup = 'SNR_Vnot0'
# method = 'main_v5_v4_multifreqtrun'
# default_setting = 'pr=0.4,n=30,beta=1.3'
# 
# for (freq_trun in c(9,7,5,3,1)) {
#   
#   for (. in 1:split) {
#     ### N_node
#     for (i in 1:length(N_node_persubj_list)) {
#       N_node = N_node_persubj_list[[i]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5(SEED = SEED, 
#                          N_node_vec = rep(N_node,1),
#                          conn_prob_mean = 0.4, 
#                          conn_patt_sep = 1.3,
#                          time_shift_mean_vec = rep(20,N_clus),
#                          t_vec = seq(0,200,length.out=200),
#                          freq_trun=freq_trun,
#                          N_clus_min = 1, N_clus_max = 5),
#                  error = function(x) print(SEED))
#       }
#       param_name = "n"
#       param_value = N_node
#       folder_path = paste0(top_level_folder, '/', setup, '/', method, '/',
#                            default_setting, '/', param_name, '/', param_value,
#                            '/', 'freqtrun','/',freq_trun)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     
#     # ### beta
#     # for (i in 1:length(conn_patt_sep_list)) {
#     #   conn_patt_sep = conn_patt_sep_list[[i]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5(SEED = SEED, 
#     #                      N_node_vec = rep(30,1),
#     #                      conn_prob_mean = 0.4, 
#     #                      conn_patt_sep = conn_patt_sep,
#     #                      time_shift_mean_vec = rep(20,N_clus),
#     #                      t_vec = seq(0,200,length.out=200),
#     #                      freq_trun=freq_trun,
#     #                      N_clus_min = 1, N_clus_max = 5),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "beta"
#     #   param_value = conn_patt_sep
#     #   folder_path = paste0(top_level_folder, '/', setup, '/', method, '/',
#     #                        default_setting, '/', param_name, '/', param_value,
#     #                        '/', 'freqtrun','/',freq_trun)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     #   
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
### SNR, V==0 -----
# ### Various signal-to-noise ratio: n, beta, clus_size_vec, V==0, alpha, p
# ### main_v3.2.1 (cdf) / main_v3.4.1 (pdf)
# ### Setup: n
# N_node_persubj_list = list(30,36,42,48,54,60,66,72,78,84,90)
# ### beta
# conn_patt_sep_list = list(1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8)
# 
# 
# ### pdf ###
# top_level_folder = "../Results/Rdata"
# setup = 'SNR_Vis0'
# method = 'our_v3.4.1'
# default_setting = 'pr=0.4,n=30,beta=1.3'
# 
# for (. in 1:split) {
#   ### N_node
#   for (i in 1:length(N_node_persubj_list)) {
#     N_node = N_node_persubj_list[[i]]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.4.1(SEED = SEED, N_node_vec = rep(N_node,1),
#                            conn_prob_mean = 0.4, conn_patt_sep = 1.3,
#                            t_vec = seq(0,200,length.out=200),
#                            time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     param_name = "n"
#     param_value = N_node
#     folder_path = paste0(top_level_folder, '/', setup, '/', method, '/', 
#                          default_setting, '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#   }
#   
#   
#   ### beta
#   for (i in 1:length(conn_patt_sep_list)) {
#     conn_patt_sep = conn_patt_sep_list[[i]]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.4.1(SEED = SEED, conn_patt_sep = conn_patt_sep,
#                            N_node_vec = rep(30,1),
#                            conn_prob_mean = 0.4,
#                            t_vec = seq(0,200,length.out=200),
#                            time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     param_name = "beta"
#     param_value = conn_patt_sep
#     folder_path = paste0(top_level_folder, '/', setup, '/', method, '/', 
#                          default_setting, '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#   }
#   
#   
# 
# }
# 
# 
# ### cdf ###
# top_level_folder = "../Results/Rdata"
# setup = 'SNR_Vis0'
# method = 'our_v3.2.1'
# default_setting = 'pr=0.4,n=30,beta=1.3'
# 
# for (. in 1:split) {
#   ### n
#   for (i in 1:length(N_node_persubj_list)) {
#     N_node = N_node_persubj_list[[i]]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.2.1(SEED = SEED, N_node_vec = rep(N_node,1),
#                            conn_prob_mean = 0.4, conn_patt_sep = 1.3,
#                            t_vec = seq(0,200,length.out=200),
#                            time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     param_name = "n"
#     param_value = N_node
#     folder_path = paste0(top_level_folder, '/', setup, '/', method, '/', 
#                          default_setting, '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#   }
#   
#   
#   ### beta
#   for (i in 1:length(conn_patt_sep_list)) {
#     conn_patt_sep = conn_patt_sep_list[[i]]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.2.1(SEED = SEED, conn_patt_sep = conn_patt_sep,
#                            N_node_vec = rep(30,1),
#                            conn_prob_mean = 0.4,
#                            t_vec = seq(0,200,length.out=200),
#                            time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     param_name = "beta"
#     param_value = conn_patt_sep
#     folder_path = paste0(top_level_folder, '/', setup, '/', method, '/', 
#                          default_setting, '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#   }
#   
#   
#   
# }
# 








### Multi-subj (Case 1-3) ----
# now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
# folder_path_1 = paste0("../Results/Rdata/", 'Multi_subj_case1_v3.1.1/n=30,beta=1.5,p=0.7/', now_sim)
# folder_path_2 = paste0("../Results/Rdata/", 'Multi_subj_case2_v3.1.1/n=30,beta=1.5,p=0.7/', now_sim)
# folder_path_3 = paste0("../Results/Rdata/", 'Multi_subj_case3_v3.1.1/n=30,beta=1.5,p=0.7/', now_sim)
# folder_path_4 = paste0("../Results/Rdata/", 'Multi_subj_case4_v3.3.1/n=30,beta=1.5,p=0.7/', now_sim)
# dir.create(path = folder_path_1, recursive = TRUE, showWarnings = FALSE)
# dir.create(path = folder_path_2, recursive = TRUE, showWarnings = FALSE)
# dir.create(path = folder_path_3, recursive = TRUE, showWarnings = FALSE)
# dir.create(path = folder_path_4, recursive = TRUE, showWarnings = FALSE)
# for (. in 1:split) {
#   
#   ### Case 1: change N_subj
#   results1 = list()
#   for (i in 1:length(N_subj_list)) {
#     N_subj = N_subj_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.1.1(SEED=SEED, N_subj=N_subj, N_node_vec = rep(30,N_subj),
#                            conn_prob_mean = 0.7, conn_patt_sep = 1.5),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("N_subj_",N_subj)]] = tmp
#   }
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path_1,'/','Case1','_', 'N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))
#   rm(results1)
#   
#   
#   ### Case 2: unbalanced network size
#   results1 = list()
#   for (i in 1:length(N_subj_list)) {
#     N_subj = N_subj_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       N_node_vec = c(rmultinom(1,9*N_subj,c(1/2,rep(1/2/(N_subj-1),N_subj-1)))*3 + 3)
#       # N_node_vec = c(18*N_subj+12,rep(12,N_subj-1))
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.1.1(SEED=SEED, N_subj=N_subj, N_node_vec = N_node_vec,
#                            conn_prob_mean = 0.7, conn_patt_sep = 1.5),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("N_subj_",N_subj)]] = tmp
#   }
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path_2,'/','Case2','_', 'N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))
#   rm(results1)
#   
#   
#   ### Case 3: unbalanced cluster size
#   results1 = list()
#   for (i in 1:length(N_subj_list)) {
#     N_subj = N_subj_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       clus_size_mat = t(replicate(N_subj, rmultinom(1,27,sample(c(1/2,1/4,1/4)))+1, simplify = 'matrix'))
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.1.1(SEED=SEED,
#                            N_subj=N_subj, N_node_vec = rep(30,N_subj),
#                            clus_size_mat = clus_size_mat,
#                            conn_prob_mean = 0.7, conn_patt_sep = 1.5),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("N_subj_",N_subj)]] = tmp
#   }
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path_3,'/','Case3','_', 'N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))
#   rm(results1)
#   
#   
#   ### Case 4: baseline --- F is given as F_true
#   results1 = list()
#   for (i in 1:length(N_subj_list)) {
#     N_subj = N_subj_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.3.1(SEED=SEED,
#                            N_subj=N_subj, N_node_vec = rep(30,N_subj),
#                            conn_prob_mean = 0.7, conn_patt_sep = 1.5),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("N_subj_",N_subj)]] = tmp
#   }
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path_4,'/','Case4','_', 'N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))
#   rm(results1)
#   
# }
# 




### Jitter_radius -----
# 
# now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
# folder_path_adp = paste0("../Results/Rdata/", 'Jitter_radius/','Adaptive_order/', now_sim)
# folder_path_fix = paste0("../Results/Rdata/", 'Jitter_radius/','Fixed_order/', now_sim)
# dir.create(path = folder_path_adp, recursive = TRUE, showWarnings = FALSE)
# dir.create(path = folder_path_fix, recursive = TRUE, showWarnings = FALSE)
# 
# for (. in 1:split) {
#   
#   results1 = list()
#   for (i in 1:length(jitter_time_rad_list)) {
#     jitter_time_rad = jitter_time_rad_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v3.1(SEED=SEED, jitter_time_rad = jitter_time_rad,
#                          conn_prob_mean = 0.5, conn_patt_sep = 1.2),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("jitter_time_rad_",jitter_time_rad)]] = tmp
#   }
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path_adp,'/','Jitter_radius_','N_trial', N_trial, '_Adaptive', '_', now_trial, '.Rdata'))
#   rm(results1)
#   
#   
#   results1 = list()
#   for (i in 1:length(jitter_time_rad_list)) {
#     jitter_time_rad = jitter_time_rad_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v4.1(SEED=SEED, jitter_time_rad = jitter_time_rad,
#                          conn_prob_mean = 0.5, conn_patt_sep = 1.2),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("jitter_time_rad_",jitter_time_rad)]] = tmp
#   }
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path_fix,'/','Jitter_radius_','N_trial', N_trial, '_Adaptive', '_', now_trial, '.Rdata'))
#   rm(results1)
# }
# 
# 





### Unbalanced clusters, single subj ----
# ##### unbalancedness
# n_star_list = list(20,24,28,32,36,40,44,48)
# 
# results1 = results2 = results3 = results4 = results5 = results6 = list()
# 
# now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
# folder_path = paste0("../Results/Rdata/", 'Unbalanced_clusters_our_v3.2/p=0.5,n=60,beta=1.5,V=0/', now_sim)
# dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
# for (. in 1:split) {
# 
#   ### (n_star,(30-n_star)/2,(30-n_star)/2)
#   for (i in 1:length(n_star_list)) {
#     n_star = n_star_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       clus_size_vec = (c(n_star,(60-n_star)/2,(60-n_star)/2))
#       tryCatch(main_v3.2(SEED = SEED, N_node_vec = rep(60,1),
#                          clus_size_mat = matrix(clus_size_vec,
#                                                  nrow=1,
#                                                  ncol=3,
#                                                  byrow = TRUE),
#                          conn_prob_mean = 0.5, conn_patt_sep = 1.5,
#                          time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("n_star_",n_star)]] = tmp
#   }
# 
# 
#   ### ((30-n_star)/2,n_star,(30-n_star)/2)
#   for (i in 1:length(n_star_list)) {
#     n_star = n_star_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       clus_size_vec = (c((60-n_star)/2,n_star,(60-n_star)/2))
#       tryCatch(main_v3.2(SEED = SEED, N_node_vec = rep(60,1),
#                          clus_size_mat = matrix(clus_size_vec,
#                                                 nrow=1,
#                                                 ncol=3,
#                                                 byrow = TRUE),
#                          conn_prob_mean = 0.5, conn_patt_sep = 1.5,
#                          time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     results2[[paste0("n_star_",n_star)]] = tmp
#   }
# 
# 
#   ### ((30-n_star)/2,(30-n_star)/2,n_star)
#   for (i in 1:length(n_star_list)) {
#     n_star = n_star_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       clus_size_vec = (c((60-n_star)/2,(60-n_star)/2,n_star))
#       tryCatch(main_v3.2(SEED = SEED, N_node_vec = rep(60,1),
#                          clus_size_mat = matrix(clus_size_vec,
#                                                 nrow=1,
#                                                 ncol=3,
#                                                 byrow = TRUE),
#                          conn_prob_mean = 0.5, conn_patt_sep = 1.5,
#                          time_shift_mean_vec = rep(0,N_clus)),
#                error = function(x) print(SEED))
#     }
#     results3[[paste0("n_star_",n_star)]] = tmp
#   }
# 
# 
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path, '/','N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))
# 
# }
# 





