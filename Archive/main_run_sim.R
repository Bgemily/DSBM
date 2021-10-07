#!/usr/bin/env Rscript

# Case 1 ---------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)



library("optparse")

option_list = list(
  make_option(c("-n", "--NSim"), type="integer", default=1, 
              help="number of repeated trials"),
  make_option(c("-c", "--case"), type="integer", default=2, 
              help="choose simulation setting"),
  make_option(c("-k", "--N_clus"), type="integer", default=3, 
              help="number of clusters"),
  make_option("--tau", type="integer", default=40),
  make_option("--split", type="integer", default=10)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


NSim = opt$NSim
case = opt$case
N_clus = opt$N_clus
tau_std = opt$tau

Ncores = 20

# step_size = 0.02



library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)


results1 = list()


######### V1
# results1 <- foreach(i = 1:NSim) %dopar% {
#   SEED = SEED_vec[i]
#   main(case=1, SEED=SEED, N_clus=2, N_overclus=3, MaxIter = 2)
# }
######## V2
tau_max_vec = c(20,40,60,80,100,120,140,160,180)
# tau_std = 40
conn_prob_vec = c(0.7,0.5,0.4,0.3,0.2,0.19,0.18,0.17,0.16,0.15,0.14,0.13)
conn_prob_std = 0.3
SEED = 831
beta_vec = c(1.6, 1.5, 1.4, 1.3, 1.2, 1.1)
beta_std = 1.3
alpha_vec = c(0.5,1,2,3,4,6,8)
alpha_std = 1
total_time = 200
clus_size_vec_list = list(c(30,30,30), c(32,32,26), c(34,34,22), c(36,36,18), 
                          c(38,38,14), c(40,40,10), c(42,42,6))
# conv_threshold_vec = c(1e-1,1e-2,1e-3,5e-4,2e-4,1e-4,1e-5,1e-6,1e-7)
kmeans = FALSE # if FALSE, use k-medoids

seed = 9012 
set.seed(seed)

NSim_total = NSim
split = opt$split
NSim = NSim_total/split
for (. in 1:split) {
  
  for (i in 1:length(tau_max_vec)) {
    tmp <- foreach(j = 1:NSim) %dopar% {
      # SEED = SEED+10
      tryCatch(main(case=case, SEED=NULL, N_clus=N_clus, MaxIter = 10,
                    tau_max = tau_max_vec[i], conn_prob = conn_prob_std,
                    beta=beta_std, alpha=alpha_std,
                    tau_struc = max,standardize = T, total_time = total_time),
               error = function(x) print(SEED))
    }
    results1[[paste0("tau_max_",tau_max_vec[i])]] = tmp
  }



  for (i in 1:length(conn_prob_vec)) {
    tmp <- foreach(j = 1:NSim) %dopar% {
      # SEED = SEED+10
      tryCatch(main(case=case, SEED=NULL, N_clus=N_clus, MaxIter = 10,
                    tau_max = tau_std, conn_prob = conn_prob_vec[i],
                    beta=beta_std, alpha=alpha_std,
                    tau_struc = max,standardize = T, total_time = total_time),
               error = function(x) print(SEED))
    }
    results1[[paste0("conn_prob_",conn_prob_vec[i])]] = tmp
  }



  for (i in 1:length(beta_vec)) {
    tmp <- foreach(j = 1:NSim) %dopar% {
      # SEED = SEED+10
      tryCatch(main(case=case, SEED=NULL, N_clus=N_clus, MaxIter = 10,
                    tau_max = tau_std, conn_prob = conn_prob_std,
                    beta=beta_vec[i], alpha=alpha_std,
                    tau_struc = max,standardize = T, total_time = total_time),
               error = function(x) print(SEED))
    }
    results1[[paste0("beta_",beta_vec[i])]] = tmp
  }



  for (i in 1:length(alpha_vec)) {
    tmp <- foreach(j = 1:NSim) %dopar% {
      # SEED = SEED+10
      tryCatch(main(case=case, SEED=NULL, N_clus=N_clus, MaxIter = 10,
                    tau_max = tau_std, conn_prob = conn_prob_std,
                    beta=beta_std, alpha=alpha_vec[i],
                    tau_struc = max,standardize = T, total_time = total_time),
               error = function(x) print(SEED))
    }
    results1[[paste0("alpha_",alpha_vec[i])]] = tmp
  }
  
  for (i in 1:length(clus_size_vec_list)) {
    tmp <- foreach(j = 1:NSim) %dopar% {
      # SEED = SEED+10
      tryCatch(main(case=case, SEED=NULL, N_clus=N_clus, MaxIter = 10,
                    tau_max = tau_std, conn_prob = conn_prob_std,
                    beta=beta_std, alpha=alpha_std,
                    tau_struc = max,standardize = T, total_time = total_time,
                    clus_size_vec = clus_size_vec_list[[i]]),
               error = function(x) print(SEED))
    }
    results1[[paste0("clus_size_",clus_size_vec_list[i])]] = tmp
  }
  
  
  
  # for (i in 1:length(conv_threshold_vec)) {
  #   tmp <- foreach(j = 1:NSim) %dopar% {
  #     # SEED = SEED+10
  #     tryCatch(main(case=case, SEED=NULL, N_clus=N_clus, MaxIter = 100,
  #                   tau_max = tau_std, conn_prob = conn_prob_std,
  #                   beta=beta_std, alpha=alpha_std,
  #                   tau_struc = max,standardize = T, total_time = total_time,
  #                   conv_thres = conv_threshold_vec[i]
  #                   ),
  #              error = function(x) print(SEED))
  #   }
  #   results1[[paste0("conv_thres_",conv_threshold_vec[i])]] = tmp
  # }
  
  #######################
  
  
  now = format(Sys.time(), "%Y%m%d_%H%M%S")
  save.image(paste0('case',case,'_NSim', NSim, '_our', '_', now, '.Rdata'))
  
  
}




###### debug
# tau_max = 20
# conn_prob = 1
# SEED_vec = seq(189,110765,length.out=100)
# i=3
# SEED = SEED_vec[i]
# beta=1
# 
# main(case=case, SEED=SEED, N_clus=N_clus, MaxIter = 10, tau_max = tau_max, 
#      conn_prob = conn_prob, beta=beta, tau_struc = max)->r
# r$clus_result$clusters_history
# true_pdf_array = fun2pdfarray(true_pdf_fun_list = r$network$true_pdf_fun_list,tau_mat = r$network$tau_mat,
#                               membership_true = r$network$membership_true,t_vec = r$network$t_vec)
# plot_pdf_array(r$clus_result$center_pdf_array, true_pdf_array)
# 
