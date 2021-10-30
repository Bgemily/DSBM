#!/usr/bin/env Rscript

### Save results with structure: setup/method/default_setting/param_names/param_values/time_stamp.Rdata
### For example, setup == "SNR_Vis0", "SNR_Vnot0", "Unbalanced_clusters", "Multi_subj"
### method == "ppsbm"
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
  make_option(c("-n", "--N_trial"), type="integer", default=20, 
              help="number of repeated trials"),
  make_option(c("s","--split"), type="integer", default=4)
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
const = 40
max_iter = 3 ### number of iterations when updating time shifts


# Run simulations ---------------------------------------------------------

### N_clus_est, V==0 -----

### Parameters' possible values: 
### n
N_node_persubj_list = list(30,42,54,66,78,90)
# N_node_persubj_list = list(90)
### beta
conn_patt_sep_list = list(1.3,1.4,1.5,1.6,1.7,1.8)
# conn_patt_sep_list = list(1.8)


### Method: apply_ppsbm_ICL ###
top_level_folder = "../Results/Rdata"
setup = 'SNR_Vis0'
method = 'apply_ppsbm_ICL'

default_setting = 'pr=0.4,n=30,beta=1.8'

for (. in 1:split) {
  ### N_node
  for (i in 1:length(N_node_persubj_list)) {
    N_node = N_node_persubj_list[[i]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(apply_ppsbm_ICL(SEED = SEED, 
                               N_node_vec = rep(N_node,1),
                               conn_prob_mean = 0.4, 
                               conn_patt_sep = 1.8,
                               time_shift_mean_vec = rep(0,3),
                               t_vec = seq(0,200,length.out=200),
                               Qmin = 1, Qmax = 5),
               error = function(x) print(SEED))
    }
    param_name = "n"
    param_value = N_node
    folder_path = paste0(top_level_folder, '/', setup, '/', method, '/',
                         default_setting, '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  
  
}

default_setting = 'pr=0.4,n=30,beta=1.3'

for (. in 1:split) {
  ### N_node
  for (i in 1:length(N_node_persubj_list)) {
    N_node = N_node_persubj_list[[i]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(apply_ppsbm_ICL(SEED = SEED, 
                               N_node_vec = rep(N_node,1),
                               conn_prob_mean = 0.4, 
                               conn_patt_sep = 1.3,
                               time_shift_mean_vec = rep(0,3),
                               t_vec = seq(0,200,length.out=200),
                               Qmin = 1, Qmax = 5),
               error = function(x) print(SEED))
    }
    param_name = "n"
    param_value = N_node
    folder_path = paste0(top_level_folder, '/', setup, '/', method, '/',
                         default_setting, '/', param_name, '/', param_value)
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
      tryCatch(apply_ppsbm_ICL(SEED = SEED, 
                               N_node_vec = rep(30,1),
                               conn_prob_mean = 0.4, 
                               conn_patt_sep = conn_patt_sep,
                               time_shift_mean_vec = rep(0,3),
                               t_vec = seq(0,200,length.out=200),
                               Qmin = 1, Qmax = 5),
               error = function(x) print(SEED))
    }
    param_name = "beta"
    param_value = conn_patt_sep
    folder_path = paste0(top_level_folder, '/', setup, '/', method, '/',
                         default_setting, '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }


}




### TODO: modify the following saving format of result.
# ### Various signal-to-noise ratio: n, beta, clus_size_vec, V!=0, alpha, p ----
# ### Setup: n
# N_node_persubj_list = list(30,36,42,48,54,60,66,72,78,84,90)
# ### beta
# conn_patt_sep_list = list(1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0)
# # ### unbalancedness
# # n_star_list = list(22,18,14,10)
# ### V
# time_shift_mean_vec_list = list(rep(15,N_clus), rep(17.5,N_clus), 
#                                 rep(20,N_clus), rep(22.5,N_clus),
#                                 rep(25,N_clus), rep(27.5,N_clus), 
#                                 rep(30,N_clus), rep(32.5,N_clus),
#                                 rep(35,N_clus), rep(37.5,N_clus), rep(40,N_clus))
# ### ### ###
# results1 = results2 = results3 = results4 = results5 = results6 = list()
# 
# now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
# folder_path = paste0("../Results/Rdata/", 'SNR_Vneq0_ppsbm_pdf_corrected_fmse/p=0.7,n=30,beta=1.5,V=80/', now_sim)
# dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
# for (. in 1:split) {
#   
#   ### n
#   for (i in 1:length(N_node_persubj_list)) {
#     N_node = N_node_persubj_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(N_node,1),
#                               conn_prob_mean = 0.7, conn_patt_sep = 1.5,
#                               t_vec = seq(0,200,length.out=200),
#                               time_shift_mean_vec = rep(40,N_clus)),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("N_node_",N_node)]] = tmp
#   }
#   
#   ### beta
#   for (i in 1:length(conn_patt_sep_list)) {
#     conn_patt_sep = conn_patt_sep_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(apply_ppsbm_v2(SEED = SEED, conn_patt_sep = conn_patt_sep,
#                               N_node_vec = rep(30,1),
#                               conn_prob_mean = 0.7,
#                               t_vec = seq(0,200,length.out=200),
#                               time_shift_mean_vec = rep(40,N_clus)),
#                error = function(x) print(SEED))
#     }
#     results2[[paste0("conn_patt_sep_",conn_patt_sep)]] = tmp
#   }
#   
#   ### n_star: clus_size_vec
#   # for (i in 1:length(n_star_list)) {
#   #   n_star = n_star_list[[i]]
#   #   tmp <- foreach(j = 1:N_trial) %dopar% {
#   #     SEED = sample(1:1e7,1)
#   #     clus_size_vec = sample(c(n_star,(60-n_star)/2,(60-n_star)/2))
#   #     tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(60,1),
#   #                             clus_size_mat = matrix(clus_size_vec,
#   #                                                    nrow=1,
#   #                                                    ncol=3,
#   #                                                    byrow = TRUE),
#   #                             conn_prob_mean = 1, conn_patt_sep = 1.5,
#   #                             time_shift_mean_vec = rep(40,N_clus)),
#   #              error = function(x) print(SEED))
#   #   }
#   #   results3[[paste0("n_star_",n_star)]] = tmp
#   # }
#   # 
#   ### V: time_shift_mean_vec_list
#   for (i in 1:length(time_shift_mean_vec_list)) {
#     time_shift_mean_vec = time_shift_mean_vec_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(apply_ppsbm_v2(SEED = SEED, time_shift_mean_vec = time_shift_mean_vec,
#                               N_node_vec = rep(30,1),
#                               t_vec = seq(0,200,length.out=200),
#                               conn_prob_mean = 0.7, conn_patt_sep = 1.5),
#                error = function(x) print(SEED))
#     }
#     results4[[paste0("time_shift_mean_vec_",(time_shift_mean_vec[1]))]] = tmp
#   }
#   
#   ### alpha: conn_patt_var_list
#   # for (i in 1:length(conn_patt_var_list)) {
#   #   conn_patt_var = conn_patt_var_list[[i]]
#   #   tmp <- foreach(j = 1:N_trial) %dopar% {
#   #     SEED = sample(1:1e7,1)
#   #     tryCatch(apply_ppsbm_v2(SEED = SEED, conn_patt_var = conn_patt_var,
#   #                             conn_prob_mean = 0.7),
#   #              error = function(x) print(SEED))
#   #   }
#   #   results5[[paste0("conn_patt_var_",conn_patt_var)]] = tmp
#   # }
#   
#   
#   ### p: conn_prob_mean_list
#   # for (i in 1:length(conn_prob_mean_list)) {
#   #   conn_prob_mean = conn_prob_mean_list[[i]]
#   #   tmp <- foreach(j = 1:N_trial) %dopar% {
#   #     SEED = sample(1:1e7,1)
#   #     tryCatch(apply_ppsbm_v2(SEED = SEED, conn_prob_mean = conn_prob_mean),
#   #              error = function(x) print(SEED))
#   #   }
#   #   results6[[paste0("conn_prob_mean_",conn_prob_mean)]] = tmp
#   # }
#   
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path, '/','N_trial', N_trial, '_ppsbm', '_', now_trial, '.Rdata'))
#   
# }
# 
# 
# 
# 
# ### Various signal-to-noise ratio: n, beta, clus_size_vec, V==0, alpha, p ----
# ### Setup: n
# N_node_persubj_list = list(30,36,42,48,54,60,66,72,78,84,90)
# ### beta
# conn_patt_sep_list = list(1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8)
# # ### unbalancedness
# # n_star_list = list(22,18,14,10)
# ### ### ###
# results1 = results2 = results3 = results4 = results5 = results6 = list()
# 
# now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
# folder_path = paste0("../Results/Rdata/", 'SNR_Vis0_ppsbm_pdf_corrected_fmse/','p=0.4,n=30,beta=1.3/', now_sim)
# dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
# for (. in 1:split) {
#   
#   ### n
#   for (i in 1:length(N_node_persubj_list)) {
#     N_node = N_node_persubj_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(N_node,1),
#                               time_shift_mean_vec = rep(0,N_clus),
#                               t_vec = seq(0,200,length.out=200),
#                               conn_prob_mean = 0.4, conn_patt_sep = 1.3),
#                error = function(x) print(SEED))
#     }
#     results1[[paste0("N_node_",N_node)]] = tmp
#   }
#   
#   ### beta
#   for (i in 1:length(conn_patt_sep_list)) {
#     conn_patt_sep = conn_patt_sep_list[[i]]
#     tmp <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(apply_ppsbm_v2(SEED = SEED, conn_patt_sep = conn_patt_sep,
#                               N_node_vec = rep(30,1),
#                               time_shift_mean_vec = rep(0,N_clus),
#                               t_vec = seq(0,200,length.out=200),
#                               conn_prob_mean = 0.4),
#                error = function(x) print(SEED))
#     }
#     results2[[paste0("conn_patt_sep_",conn_patt_sep)]] = tmp
#   }
#   
#   ### n_star: clus_size_vec
#   # for (i in 1:length(n_star_list)) {
#   #   n_star = n_star_list[[i]]
#   #   tmp <- foreach(j = 1:N_trial) %dopar% {
#   #     SEED = sample(1:1e7,1)
#   #     clus_size_vec = sample(c(n_star,(60-n_star)/2,(60-n_star)/2))
#   #     tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(60,1),
#   #                             clus_size_mat = matrix(clus_size_vec,
#   #                                                    nrow=1,
#   #                                                    ncol=3,
#   #                                                    byrow = TRUE),
#   #                             time_shift_mean_vec = rep(0,N_clus),
#   #                             conn_prob_mean = 0.3, conn_patt_sep = 1.25),
#   #              error = function(x) print(SEED))
#   #   }
#   #   results3[[paste0("n_star_",n_star)]] = tmp
#   # }
#   # 
#   
#   ### alpha: conn_patt_var_list
#   # for (i in 1:length(conn_patt_var_list)) {
#   #   conn_patt_var = conn_patt_var_list[[i]]
#   #   tmp <- foreach(j = 1:N_trial) %dopar% {
#   #     SEED = sample(1:1e7,1)
#   #     tryCatch(apply_ppsbm_v2(SEED = SEED, conn_patt_var = conn_patt_var,
#   #                             time_shift_mean_vec = rep(0,N_clus),
#   #                             conn_prob_mean = 0.5),
#   #              error = function(x) print(SEED))
#   #   }
#   #   results5[[paste0("conn_patt_var_",conn_patt_var)]] = tmp
#   # }
#   
#   
#   ### p: conn_prob_mean_list
#   # for (i in 1:length(conn_prob_mean_list)) {
#   #   conn_prob_mean = conn_prob_mean_list[[i]]
#   #   tmp <- foreach(j = 1:N_trial) %dopar% {
#   #     SEED = sample(1:1e7,1)
#   #     tryCatch(apply_ppsbm_v2(SEED = SEED, conn_prob_mean = conn_prob_mean,
#   #                             time_shift_mean_vec = rep(0,N_clus),
#   #                             conn_prob_mean = 0.5),
#   #              error = function(x) print(SEED))
#   #   }
#   #   results6[[paste0("conn_prob_mean_",conn_prob_mean)]] = tmp
#   # }
#   
#   
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save.image(paste0(folder_path, '/','N_trial', N_trial, '_ppsbm', '_', now_trial, '.Rdata'))
#   
# }
# 
# 
# # 
# # ### Unbalanced clusters, single subj
# # ##### unbalancedness
# # n_star_list = list(20,24,28,32,36,40,44,48)
# # 
# # results1 = results2 = results3 = results4 = results5 = results6 = list()
# # 
# # now_sim = format(Sys.time(), "%Y%m%d_%H%M%S")
# # folder_path = paste0("../Results/Rdata/", 'Unbalanced_clusters_ppsbm/p=0.5,n=60,beta=1.5,V=0/', now_sim)
# # dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
# # for (. in 1:split) {
# #   
# #   ### (n_star,(60-n_star)/2,(60-n_star)/2)
# #   for (i in 1:length(n_star_list)) {
# #     n_star = n_star_list[[i]]
# #     tmp <- foreach(j = 1:N_trial) %dopar% {
# #       SEED = sample(1:1e7,1)
# #       clus_size_vec = (c(n_star,(60-n_star)/2,(60-n_star)/2))
# #       tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(60,1),
# #                          clus_size_mat = matrix(clus_size_vec,
# #                                                 nrow=1,
# #                                                 ncol=3,
# #                                                 byrow = TRUE),
# #                          conn_prob_mean = 0.5, conn_patt_sep = 1.5,
# #                          time_shift_mean_vec = rep(0,N_clus)),
# #                error = function(x) print(SEED))
# #     }
# #     results1[[paste0("n_star_",n_star)]] = tmp
# #   }
# #   
# #   
# #   ### ((60-n_star)/2,n_star,(60-n_star)/2)
# #   for (i in 1:length(n_star_list)) {
# #     n_star = n_star_list[[i]]
# #     tmp <- foreach(j = 1:N_trial) %dopar% {
# #       SEED = sample(1:1e7,1)
# #       clus_size_vec = (c((60-n_star)/2,n_star,(60-n_star)/2))
# #       tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(60,1),
# #                          clus_size_mat = matrix(clus_size_vec,
# #                                                 nrow=1,
# #                                                 ncol=3,
# #                                                 byrow = TRUE),
# #                          conn_prob_mean = 0.5, conn_patt_sep = 1.5,
# #                          time_shift_mean_vec = rep(0,N_clus)),
# #                error = function(x) print(SEED))
# #     }
# #     results2[[paste0("n_star_",n_star)]] = tmp
# #   }
# #   
# #   
# #   ### ((60-n_star)/2,(60-n_star)/2,n_star)
# #   for (i in 1:length(n_star_list)) {
# #     n_star = n_star_list[[i]]
# #     tmp <- foreach(j = 1:N_trial) %dopar% {
# #       SEED = sample(1:1e7,1)
# #       clus_size_vec = (c((60-n_star)/2,(60-n_star)/2,n_star))
# #       tryCatch(apply_ppsbm_v2(SEED = SEED, N_node_vec = rep(60,1),
# #                          clus_size_mat = matrix(clus_size_vec,
# #                                                 nrow=1,
# #                                                 ncol=3,
# #                                                 byrow = TRUE),
# #                          conn_prob_mean = 0.5, conn_patt_sep = 1.5,
# #                          time_shift_mean_vec = rep(0,N_clus)),
# #                error = function(x) print(SEED))
# #     }
# #     results3[[paste0("n_star_",n_star)]] = tmp
# #   }
# #   
# #   
# #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
# #   save.image(paste0(folder_path, '/','N_trial', N_trial, '_our', '_', now_trial, '.Rdata'))
# #   
# # }
# 
# 
# 



