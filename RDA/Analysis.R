
# Default working directory: Sources/

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)
library(scales)
library(combinat)
library(cluster)
library(mclust)
library(reshape2)
library(ppsbm)
library(foreach)
library(doParallel)


# Parallel computing setup ------------------------------

N_cores = 20
registerDoParallel(cores=N_cores)


# Get network information for the L/R side ------------------------------------

data_folder = "../Processed_FunctionalData/"
path_vec = list.files(data_folder, full.names = TRUE)
file_vec = list.files(data_folder, full.names = FALSE)
edge_time_mat_list = vector(mode = "list", length = length(path_vec)*2)
for(m in 1:length(path_vec)){ 
  path = path_vec[m]
  
  # ### Read information from data
  # edge_time_mat_L = as.matrix(read.csv(paste(path, '/EdgeTime_L.csv', sep='')))
  # edge_time_mat_L = edge_time_mat_L[,-1]
  # edge_time_mat_R = as.matrix(read.csv(paste(path, '/EdgeTime_R.csv', sep='')))
  # edge_time_mat_R = edge_time_mat_R[,-1]
  # 
  # edge_time_mat_list[c(2*m-1,2*m)] = list(edge_time_mat_L, edge_time_mat_R)
  
  ### Remove neurons with large coordinate along DV axis
  edge_time_mat = as.matrix(read.csv(paste(path, '/EdgeTime.csv', sep='')))
  edge_time_mat = edge_time_mat[,-1]
  avai.inds = as.matrix(read.csv(paste(path,'/AvaiNeurons.csv',sep='')))
  avai.inds = avai.inds[,-1];
  locs.all = as.matrix(read.csv(paste(path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs_mat = locs.all[avai.inds,]
  inds_L = which(locs_mat[,2]<0)
  inds_R = which(locs_mat[,2]>0)
  locs_mat_L = locs_mat[inds_L,]
  locs_mat_R = locs_mat[inds_R,]
  edge_time_mat_L = edge_time_mat[inds_L, inds_L]
  edge_time_mat_R = edge_time_mat[inds_R, inds_R]
  avai_inds_L = which(rowSums(edge_time_mat_L<Inf)>0 & locs_mat_L[,3]<=64)
  avai_inds_R = which(rowSums(edge_time_mat_R<Inf)>0 & locs_mat_R[,3]<=64)
  edge_time_mat_L = edge_time_mat_L[avai_inds_L, avai_inds_L]
  edge_time_mat_R = edge_time_mat_R[avai_inds_R, avai_inds_R]

  edge_time_mat_list[c(2*m-1,2*m)] = list(edge_time_mat_L, edge_time_mat_R)
  
}


# Apply algorithm (our) ---------------------------------------------------------

N_clus_min = 4 # Number of clusters
N_clus_max = 4
MaxIter = 10 # Maximal iteration number
# bw = 5 # Smoothing bandwidth
conv_thres=1e-3
max_iter=10
# step_size = 0.5
N_restart = 20

max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
total_time = max_time + 10
t_vec = seq(0, total_time, 1)

N_clus_tmp = 4

### Change 4 init memberships from small clusters to rest clusters
for (subj in 4:4) {
  compl_log_lik_best = -Inf
  res_best = NULL
  # seed_restart = sample(1e5,1)
  # set.seed(seed_restart)
  fitmodel_tmp = function(ind_restart){
    ### Apply method with freq_trun and total_time_cutoff
    edge_time_mat_list_tmp = edge_time_mat_list[subj]
    res_list_tmp = list(list())
    for (ind_freq_trun_tmp in 1:1) {
      ### Get initialization
      res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
                        N_clus = N_clus_tmp,
                        t_vec = t_vec)

      clusters_list_init = res$clusters_list

      ### Add noise to initial clusters
      seed_init = sample(1e5,1)
      set.seed(seed_init)
      if (ind_restart>1) {
        # res = get_init_v5(edge_time_mat_list = edge_time_mat_list_tmp,
        #                   N_clus = N_clus_tmp,
        #                   N_restart = 1,
        #                   t_vec = t_vec)
        # clusters_list_init = res$clusters_list

        mem_init = clus2mem(clusters = clusters_list_init[[1]])
        clus_size = sapply(clusters_list_init[[1]], length)
        ind = clusters_list_init[[1]][[which.min(clus_size)]]
        ind = c(ind, sample(clusters_list_init[[1]][[which(clus_size==sort(clus_size)[2])]], 4-length(ind)))
        for (i_tmp in ind) {
          mem_init[i_tmp] = sample(setdiff(1:length(clus_size),mem_init[i_tmp]),1)
        }
        clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
      }
      
      n0_vec_list_init = res$n0_vec_list
      n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)

      # Apply algorithm
      ### Estimation z,v,f based on cdf
      time_start = Sys.time()
      res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
                            N_clus = N_clus_tmp,
                            clusters_list_init = clusters_list_init,
                            n0_vec_list_init = n0_vec_list_init,
                            n0_mat_list_init = n0_mat_list_init,
                            total_time = total_time,
                            max_iter=max_iter,
                            t_vec=t_vec,
                            freq_trun = Inf,
                            conv_thres=conv_thres,
                            MaxIter=MaxIter)
      time_end = Sys.time()
      time_estimation = time_end - time_start
      N_iteration = res$N_iteration
      res$seed_init = seed_init
      res$t_vec=t_vec
      res$clusters_list_init = clusters_list_init

      # Save results of N_clus_tmp
      res_list_tmp[[1]][[1]] = res

    }
    ### Calculate loglik
    sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
                               N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
                               N_clus_min = N_clus_tmp,
                               N_clus_max = N_clus_tmp,
                               result_list = res_list_tmp,
                               total_time = total_time)

    compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
    penalty_vec = sel_mod_res$penalty_vec

    res=res_list_tmp[[1]][[1]]
    res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec

    return(list(res=res, compl_log_lik_vec=compl_log_lik_vec,
                clusters_list_init=res$clusters_list_init))
  }
  results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
    fitmodel_tmp(ind_restart = ind_restart)
  }
  ### Update compl_log_lik_best and res_best
  compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
  clus_init_list = lapply(results_restart, '[[','clusters_list_init')

  ### Save result with freq_trun and total_time_cutoff
  method = paste0("CDF_v4_4SmallToRest",
                  "_Nrestart",N_restart,
                  "_","totaltime",total_time)
  folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
  now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
  save(results_restart,compl_log_lik_vec,clus_init_list,
       # seed_restart,
       file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))

}

# ### Change 4 init memberships from big clusters to rest clusters
# for (subj in 4:4) {
#   compl_log_lik_best = -Inf
#   res_best = NULL
#   # seed_restart = sample(1e5,1)
#   # set.seed(seed_restart)
#   fitmodel_tmp = function(ind_restart){
#     ### Apply method with freq_trun and total_time_cutoff
#     edge_time_mat_list_tmp = edge_time_mat_list[subj]
#     res_list_tmp = list(list())
#     for (ind_freq_trun_tmp in 1:1) {
#       ### Get initialization
#       res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
#                         N_clus = N_clus_tmp,
#                         t_vec = t_vec)
#       
#       clusters_list_init = res$clusters_list
#       
#       ### Add noise to initial clusters
#       seed_init = sample(1e5,1)
#       set.seed(seed_init)
#       if (ind_restart>1) {
#         # res = get_init_v5(edge_time_mat_list = edge_time_mat_list_tmp,
#         #                   N_clus = N_clus_tmp,
#         #                   N_restart = 1,
#         #                   t_vec = t_vec)
#         # clusters_list_init = res$clusters_list
#         
#         mem_init = clus2mem(clusters = clusters_list_init[[1]])
#         clus_size = sapply(clusters_list_init[[1]], length)
#         ind = sample(unlist(clusters_list_init[[1]][clus_size>4]),4)
#         mem_init[ind] = (mem_init[ind]+
#                            sample(1:(N_clus_tmp-1),length(ind),replace=TRUE)) %% N_clus_tmp
#         mem_init[mem_init==0] = N_clus_tmp
#         clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
#       }
#       
#       n0_vec_list_init = res$n0_vec_list
#       n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
#       
#       # Apply algorithm
#       ### Estimation z,v,f based on cdf
#       time_start = Sys.time()
#       res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
#                             N_clus = N_clus_tmp,
#                             clusters_list_init = clusters_list_init,
#                             n0_vec_list_init = n0_vec_list_init,
#                             n0_mat_list_init = n0_mat_list_init,
#                             total_time = total_time,
#                             max_iter=max_iter,
#                             t_vec=t_vec,
#                             freq_trun = Inf,
#                             conv_thres=conv_thres,
#                             MaxIter=MaxIter)
#       time_end = Sys.time()
#       time_estimation = time_end - time_start
#       N_iteration = res$N_iteration
#       res$seed_init = seed_init
#       res$t_vec=t_vec
#       res$clusters_list_init = clusters_list_init
#       
#       # Save results of N_clus_tmp
#       res_list_tmp[[1]][[1]] = res
#       
#     }
#     ### Calculate loglik
#     sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
#                                N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
#                                N_clus_min = N_clus_tmp,
#                                N_clus_max = N_clus_tmp,
#                                result_list = res_list_tmp,
#                                total_time = total_time)
#     
#     compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     penalty_vec = sel_mod_res$penalty_vec
#     
#     res=res_list_tmp[[1]][[1]]
#     res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     
#     return(list(res=res, compl_log_lik_vec=compl_log_lik_vec, 
#                 clusters_list_init=res$clusters_list_init))
#   }
#   results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
#     fitmodel_tmp(ind_restart = ind_restart)
#   }
#   ### Update compl_log_lik_best and res_best
#   compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
#   clus_init_list = lapply(results_restart, '[[','clusters_list_init')
#   
#   ### Save result with freq_trun and total_time_cutoff
#   method = paste0("CDF_v4_4BigToRest", 
#                   "_Nrestart",N_restart,
#                   "_","totaltime",total_time)
#   folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(results_restart,compl_log_lik_vec,clus_init_list,
#        # seed_restart, 
#        file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
#   
# }
# 
# ### Change 4 init memberships from big clusters to big clusters, i.e. not touch small clusters
# for (subj in 4:4) {
#   compl_log_lik_best = -Inf
#   res_best = NULL
#   # seed_restart = sample(1e5,1)
#   # set.seed(seed_restart)
#   fitmodel_tmp = function(ind_restart){
#     ### Apply method with freq_trun and total_time_cutoff
#     edge_time_mat_list_tmp = edge_time_mat_list[subj]
#     res_list_tmp = list(list())
#     for (ind_freq_trun_tmp in 1:1) {
#       ### Get initialization
#       res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
#                         N_clus = N_clus_tmp,
#                         t_vec = t_vec)
#       
#       clusters_list_init = res$clusters_list
#       
#       ### Add noise to initial clusters
#       seed_init = sample(1e5,1)
#       set.seed(seed_init)
#       if (ind_restart>1) {
#         mem_init = clus2mem(clusters = clusters_list_init[[1]])
#         clus_size = sapply(clusters_list_init[[1]], length)
#         ind = sample(unlist(clusters_list_init[[1]][clus_size>4]),4)
#         for (i_tmp in ind) {
#           mem_init[i_tmp] = sample(setdiff(which(clus_size>4),mem_init[i_tmp]),1)
#         }
#         clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
#       }
#       
#       n0_vec_list_init = res$n0_vec_list
#       n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
#       
#       # Apply algorithm
#       ### Estimation z,v,f based on cdf
#       time_start = Sys.time()
#       res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
#                             N_clus = N_clus_tmp,
#                             clusters_list_init = clusters_list_init,
#                             n0_vec_list_init = n0_vec_list_init,
#                             n0_mat_list_init = n0_mat_list_init,
#                             total_time = total_time,
#                             max_iter=max_iter,
#                             t_vec=t_vec,
#                             freq_trun = Inf,
#                             conv_thres=conv_thres,
#                             MaxIter=MaxIter)
#       time_end = Sys.time()
#       time_estimation = time_end - time_start
#       N_iteration = res$N_iteration
#       res$seed_init = seed_init
#       res$t_vec=t_vec
#       res$clusters_list_init = clusters_list_init
#       
#       # Save results of N_clus_tmp
#       res_list_tmp[[1]][[1]] = res
#       
#     }
#     ### Calculate loglik
#     sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
#                                N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
#                                N_clus_min = N_clus_tmp,
#                                N_clus_max = N_clus_tmp,
#                                result_list = res_list_tmp,
#                                total_time = total_time)
#     
#     compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     penalty_vec = sel_mod_res$penalty_vec
#     
#     res=res_list_tmp[[1]][[1]]
#     res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     
#     return(list(res=res, compl_log_lik_vec=compl_log_lik_vec, 
#                 clusters_list_init=res$clusters_list_init))
#   }
#   results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
#     fitmodel_tmp(ind_restart = ind_restart)
#   }
#   ### Update compl_log_lik_best and res_best
#   compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
#   clus_init_list = lapply(results_restart, '[[','clusters_list_init')
#   
#   ### Save result with freq_trun and total_time_cutoff
#   method = paste0("CDF_v4_4BigToBig", 
#                   "_Nrestart",N_restart,
#                   "_","totaltime",total_time)
#   folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(results_restart,compl_log_lik_vec,clus_init_list,
#        # seed_restart, 
#        file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
#   
# }
# 
# ### Change 8 init memberships from big clusters, i.e. not touch small clusters
# for (subj in 4:4) {
#   compl_log_lik_best = -Inf
#   res_best = NULL
#   # seed_restart = sample(1e5,1)
#   # set.seed(seed_restart)
#   fitmodel_tmp = function(ind_restart){
#     ### Apply method with freq_trun and total_time_cutoff
#     edge_time_mat_list_tmp = edge_time_mat_list[subj]
#     res_list_tmp = list(list())
#     for (ind_freq_trun_tmp in 1:1) {
#       ### Get initialization
#       res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
#                         N_clus = N_clus_tmp,
#                         t_vec = t_vec)
#       
#       clusters_list_init = res$clusters_list
#       
#       ### Add noise to initial clusters
#       seed_init = sample(1e5,1)
#       set.seed(seed_init)
#       if (ind_restart>1) {
#         mem_init = clus2mem(clusters = clusters_list_init[[1]])
#         clus_size = sapply(clusters_list_init[[1]], length)
#         ind = sample(unlist(clusters_list_init[[1]][clus_size>4]),8)
#         for (i_tmp in ind) {
#           mem_init[i_tmp] = sample(setdiff(which(clus_size>4),mem_init[i_tmp]),1)
#         }
#         while(length(unique(mem_init))<N_clus_tmp){
#           mem_init = clus2mem(clusters = clusters_list_init[[1]])
#           clus_size = sapply(clusters_list_init[[1]], length)
#           ind = sample(unlist(clusters_list_init[[1]][clus_size>4]),8)
#           for (i_tmp in ind) {
#             mem_init[i_tmp] = sample(setdiff(which(clus_size>4),mem_init[i_tmp]),1)
#           }
#         }
#         clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
#       }
#       
#       n0_vec_list_init = res$n0_vec_list
#       n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
#       
#       # Apply algorithm
#       ### Estimation z,v,f based on cdf
#       time_start = Sys.time()
#       res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
#                             N_clus = N_clus_tmp,
#                             clusters_list_init = clusters_list_init,
#                             n0_vec_list_init = n0_vec_list_init,
#                             n0_mat_list_init = n0_mat_list_init,
#                             total_time = total_time,
#                             max_iter=max_iter,
#                             t_vec=t_vec,
#                             freq_trun = Inf,
#                             conv_thres=conv_thres,
#                             MaxIter=MaxIter)
#       time_end = Sys.time()
#       time_estimation = time_end - time_start
#       N_iteration = res$N_iteration
#       res$seed_init = seed_init
#       res$t_vec=t_vec
#       res$clusters_list_init = clusters_list_init
#       
#       # Save results of N_clus_tmp
#       res_list_tmp[[1]][[1]] = res
#       
#     }
#     ### Calculate loglik
#     sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
#                                N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
#                                N_clus_min = N_clus_tmp,
#                                N_clus_max = N_clus_tmp,
#                                result_list = res_list_tmp,
#                                total_time = total_time)
#     
#     compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     penalty_vec = sel_mod_res$penalty_vec
#     
#     res=res_list_tmp[[1]][[1]]
#     res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     
#     return(list(res=res, compl_log_lik_vec=compl_log_lik_vec, 
#                 clusters_list_init=res$clusters_list_init))
#   }
#   results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
#     fitmodel_tmp(ind_restart = ind_restart)
#   }
#   ### Update compl_log_lik_best and res_best
#   compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
#   clus_init_list = lapply(results_restart, '[[','clusters_list_init')
#   
#   ### Save result with freq_trun and total_time_cutoff
#   method = paste0("CDF_v4_8BigToBig", 
#                   "_Nrestart",N_restart,
#                   "_","totaltime",total_time)
#   folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(results_restart,compl_log_lik_vec,clus_init_list,
#        # seed_restart, 
#        file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
#   
# }

# ### Pertub 4 init memberships from all clusters
# for (subj in 4:4) {
#   compl_log_lik_best = -Inf
#   res_best = NULL
#   # seed_restart = sample(1e5,1)
#   # set.seed(seed_restart)
#   fitmodel_tmp = function(ind_restart){
#     ### Apply method with freq_trun and total_time_cutoff
#     edge_time_mat_list_tmp = edge_time_mat_list[subj]
#     res_list_tmp = list(list())
#     for (ind_freq_trun_tmp in 1:1) {
#       ### Get initialization
#       res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
#                         N_clus = N_clus_tmp,
#                         t_vec = t_vec)
#       
#       clusters_list_init = res$clusters_list
#       
#       ### Add noise to initial clusters
#       seed_init = sample(1e5,1)
#       set.seed(seed_init)
#       if (ind_restart>1) {
#         mem_init = clus2mem(clusters = clusters_list_init[[1]])
#         clus_size = sapply(clusters_list_init[[1]], length)
#         ind = sample(unlist(clusters_list_init[[1]]),4)
#         pertb_mem = sample(mem_init[ind])
#         mem_init[ind] = pertb_mem
#         clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
#       }
#       
#       n0_vec_list_init = res$n0_vec_list
#       n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
# 
#       # Apply algorithm
#       ### Estimation z,v,f based on cdf
#       time_start = Sys.time()
#       res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
#                             N_clus = N_clus_tmp,
#                             clusters_list_init = clusters_list_init,
#                             n0_vec_list_init = n0_vec_list_init,
#                             n0_mat_list_init = n0_mat_list_init,
#                             total_time = total_time,
#                             max_iter=max_iter,
#                             t_vec=t_vec,
#                             freq_trun = Inf,
#                             conv_thres=conv_thres,
#                             MaxIter=MaxIter)
#       time_end = Sys.time()
#       time_estimation = time_end - time_start
#       N_iteration = res$N_iteration
#       res$seed_init = seed_init
#       res$t_vec=t_vec
#       res$clusters_list_init = clusters_list_init
#       
#       # Save results of N_clus_tmp
#       res_list_tmp[[1]][[1]] = res
#       
#     }
#     ### Calculate loglik
#     sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
#                                N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
#                                N_clus_min = N_clus_tmp,
#                                N_clus_max = N_clus_tmp,
#                                result_list = res_list_tmp,
#                                total_time = total_time)
#     
#     compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     penalty_vec = sel_mod_res$penalty_vec
#     
#     res=res_list_tmp[[1]][[1]]
#     res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     
#     return(list(res=res, compl_log_lik_vec=compl_log_lik_vec, 
#                 clusters_list_init=res$clusters_list_init))
#   }
#   results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
#     fitmodel_tmp(ind_restart = ind_restart)
#   }
#   ### Update compl_log_lik_best and res_best
#   compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
#   clus_init_list = lapply(results_restart, '[[','clusters_list_init')
#   
#   ### Save result with freq_trun and total_time_cutoff
#   method = paste0("CDF_v4_4PermAny", 
#                   "_Nrestart",N_restart,
#                   "_","totaltime",total_time)
#   folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(results_restart,compl_log_lik_vec,clus_init_list,
#        # seed_restart, 
#        file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
#   
# }


# ### Pertub 20 init memberships from all clusters
# for (subj in 4:4) {
#   compl_log_lik_best = -Inf
#   res_best = NULL
#   # seed_restart = sample(1e5,1)
#   # set.seed(seed_restart)
#   fitmodel_tmp = function(ind_restart){
#     ### Apply method with freq_trun and total_time_cutoff
#     edge_time_mat_list_tmp = edge_time_mat_list[subj]
#     res_list_tmp = list(list())
#     for (ind_freq_trun_tmp in 1:1) {
#       ### Get initialization
#       res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
#                         N_clus = N_clus_tmp,
#                         t_vec = t_vec)
#       
#       clusters_list_init = res$clusters_list
#       
#       ### Add noise to initial clusters
#       seed_init = sample(1e5,1)
#       set.seed(seed_init)
#       if (ind_restart>1) {
#         mem_init = clus2mem(clusters = clusters_list_init[[1]])
#         clus_size = sapply(clusters_list_init[[1]], length)
#         ind = sample(unlist(clusters_list_init[[1]]),20)
#         pertb_mem = sample(mem_init[ind])
#         mem_init[ind] = pertb_mem
#         clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
#       }
#       
#       n0_vec_list_init = res$n0_vec_list
#       n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
#       
#       # Apply algorithm
#       ### Estimation z,v,f based on cdf
#       time_start = Sys.time()
#       res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
#                             N_clus = N_clus_tmp,
#                             clusters_list_init = clusters_list_init,
#                             n0_vec_list_init = n0_vec_list_init,
#                             n0_mat_list_init = n0_mat_list_init,
#                             total_time = total_time,
#                             max_iter=max_iter,
#                             t_vec=t_vec,
#                             freq_trun = Inf,
#                             conv_thres=conv_thres,
#                             MaxIter=MaxIter)
#       time_end = Sys.time()
#       time_estimation = time_end - time_start
#       N_iteration = res$N_iteration
#       res$seed_init = seed_init
#       res$t_vec=t_vec
#       res$clusters_list_init = clusters_list_init
#       
#       # Save results of N_clus_tmp
#       res_list_tmp[[1]][[1]] = res
#       
#     }
#     ### Calculate loglik
#     sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
#                                N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
#                                N_clus_min = N_clus_tmp,
#                                N_clus_max = N_clus_tmp,
#                                result_list = res_list_tmp,
#                                total_time = total_time)
#     
#     compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     penalty_vec = sel_mod_res$penalty_vec
#     
#     res=res_list_tmp[[1]][[1]]
#     res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     
#     return(list(res=res, compl_log_lik_vec=compl_log_lik_vec, 
#                 clusters_list_init=res$clusters_list_init))
#   }
#   results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
#     fitmodel_tmp(ind_restart = ind_restart)
#   }
#   ### Update compl_log_lik_best and res_best
#   compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
#   clus_init_list = lapply(results_restart, '[[','clusters_list_init')
#   
#   ### Save result with freq_trun and total_time_cutoff
#   method = paste0("CDF_v4_20PermAny", 
#                   "_Nrestart",N_restart,
#                   "_","totaltime",total_time)
#   folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(results_restart,compl_log_lik_vec,clus_init_list,
#        # seed_restart, 
#        file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
#   
# }
# 
# 
# ### Pertub all init memberships from all clusters
# for (subj in 4:4) {
#   compl_log_lik_best = -Inf
#   res_best = NULL
#   # seed_restart = sample(1e5,1)
#   # set.seed(seed_restart)
#   fitmodel_tmp = function(ind_restart){
#     ### Apply method with freq_trun and total_time_cutoff
#     edge_time_mat_list_tmp = edge_time_mat_list[subj]
#     res_list_tmp = list(list())
#     for (ind_freq_trun_tmp in 1:1) {
#       ### Get initialization
#       res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp,
#                         N_clus = N_clus_tmp,
#                         t_vec = t_vec)
#       
#       clusters_list_init = res$clusters_list
#       
#       ### Add noise to initial clusters
#       seed_init = sample(1e5,1)
#       set.seed(seed_init)
#       if (ind_restart>1) {
#         mem_init = clus2mem(clusters = clusters_list_init[[1]])
#         clus_size = sapply(clusters_list_init[[1]], length)
#         mem_init = sample(mem_init)
#         clusters_list_init[[1]] = mem2clus(membership = mem_init, N_clus_min = N_clus_tmp)
#       }
#       
#       n0_vec_list_init = res$n0_vec_list
#       n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
#       
#       # Apply algorithm
#       ### Estimation z,v,f based on cdf
#       time_start = Sys.time()
#       res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list_tmp,
#                             N_clus = N_clus_tmp,
#                             clusters_list_init = clusters_list_init,
#                             n0_vec_list_init = n0_vec_list_init,
#                             n0_mat_list_init = n0_mat_list_init,
#                             total_time = total_time,
#                             max_iter=max_iter,
#                             t_vec=t_vec,
#                             freq_trun = Inf,
#                             conv_thres=conv_thres,
#                             MaxIter=MaxIter)
#       time_end = Sys.time()
#       time_estimation = time_end - time_start
#       N_iteration = res$N_iteration
#       res$seed_init = seed_init
#       res$t_vec=t_vec
#       res$clusters_list_init = clusters_list_init
#       
#       # Save results of N_clus_tmp
#       res_list_tmp[[1]][[1]] = res
#       
#     }
#     ### Calculate loglik
#     sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
#                                N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
#                                N_clus_min = N_clus_tmp,
#                                N_clus_max = N_clus_tmp,
#                                result_list = res_list_tmp,
#                                total_time = total_time)
#     
#     compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     penalty_vec = sel_mod_res$penalty_vec
#     
#     res=res_list_tmp[[1]][[1]]
#     res$compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
#     
#     return(list(res=res, compl_log_lik_vec=compl_log_lik_vec, 
#                 clusters_list_init=res$clusters_list_init))
#   }
#   results_restart <- foreach(ind_restart = 1:N_restart) %dopar% {
#     fitmodel_tmp(ind_restart = ind_restart)
#   }
#   ### Update compl_log_lik_best and res_best
#   compl_log_lik_vec = sapply(results_restart, '[[','compl_log_lik_vec')
#   clus_init_list = lapply(results_restart, '[[','clusters_list_init')
#   
#   ### Save result with freq_trun and total_time_cutoff
#   method = paste0("CDF_v4_AllPermAny", 
#                   "_Nrestart",N_restart,
#                   "_","totaltime",total_time)
#   folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(results_restart,compl_log_lik_vec,clus_init_list,
#        # seed_restart, 
#        file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
#   
# }


# # Apply algorithm (ppsbm) ---------------------------------------------------------
# 
# method = "ppsbm_saveICL"
# 
# N_clus_min = 1 # There might be errors if N_clus_min is not 1, due to bugs in ppsbm package.
# N_clus_max = 5
# max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
# total_time = max_time + 10 # Total observation length
# t_vec = seq(0, total_time, 1) # Time grid on which estimated intensities are recorded
# 
# 
# for (subj in 1:length(edge_time_mat_list)) {
#   edge_time_mat_tmp = edge_time_mat_list[[subj]]
# 
#   ### Convert data format to align with mainVEM()'s input format
#   time.seq = numeric(sum(edge_time_mat_tmp<Inf))
#   type.seq = numeric(sum(edge_time_mat_tmp<Inf))
#   current_ind = 1
#   for (i in 1:nrow(edge_time_mat_tmp)) {
#     for (j in i:ncol(edge_time_mat_tmp)) {
#       if (edge_time_mat_tmp[i,j]<Inf){
#         time.seq[current_ind] = edge_time_mat_tmp[i,j]
#         type.seq[current_ind] = convertNodePair(i, j, n = nrow(edge_time_mat_tmp), directed = FALSE)
#         current_ind = current_ind+1
#       }
#     }
#   }
#   data = list(time.seq=time.seq, type.seq=type.seq, Time=total_time)
#   Nijk = statistics(data, nrow(edge_time_mat_tmp), K=2^8, directed = FALSE)
#   data_ppsbm = list(Nijk=Nijk, Time=total_time)
# 
#   ### Apply ppsbm
#   res = mainVEM(data=data_ppsbm, n=nrow(edge_time_mat_tmp),
#                 d_part=8, Qmin=N_clus_min, Qmax=N_clus_max,
#                 directed=FALSE, method="hist")
#   if (N_clus_min == N_clus_max) {
#     res = res[[1]]
#   }
#   else {
#     sol.selec_Q = modelSelection_Q(data=data_ppsbm,n=nrow(edge_time_mat_tmp),
#                                     Qmin=N_clus_min, Qmax=N_clus_max,
#                                     directed=FALSE, sol.hist.sauv=res)
#     res = sol.selec_Q$sol.Qbest
#   }
# 
#   ### Extract estimated clusters
#   clusters_list_est = mem2clus(apply(res$tau, 2, which.max),N_clus_min=1)
# 
#   ### Extract estimated intensities
#   N_clus_est = length(clusters_list_est)
#   center_pdf_array_est = array(dim = c(N_clus_est,N_clus_est,length(t_vec)))
#   ind_qk = 1
#   for (q in 1:N_clus_est) {
#     for (k in q:N_clus_est) {
#       ### Extract estimated intensity for cluster pair (q,k)
#       intensity_tmp = exp(res$logintensities.ql[ind_qk, ])
#       ### Add breakpoints in time grid
#       rep_time = floor( (1:length(intensity_tmp))*length(t_vec)/length(intensity_tmp) ) -
#                   floor( (0:(length(intensity_tmp)-1))*length(t_vec)/length(intensity_tmp) )
#       intensity_tmp = rep(intensity_tmp, time = rep_time)
# 
#       center_pdf_array_est[q,k, ] = intensity_tmp
#       center_pdf_array_est[k,q, ] = center_pdf_array_est[q,k, ]
#       ind_qk = ind_qk + 1
#     }
#   }
# 
#   ### Save result
#   res = list(ICL_vec = sol.selec_Q$all.ICL,
#              compl_log_lik_vec = sol.selec_Q$all.compl.log.likelihood,
#              penalty_vec = sol.selec_Q$all.pen,
#              edge_time_mat = edge_time_mat_tmp,
#              clusters_list = clusters_list_est,
#              center_pdf_array = center_pdf_array_est,
#              t_vec = t_vec)
#   folder_path = paste0('../Results/Rdata/RDA/', method, '/', file_vec[(subj+1)%/%2])
#   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#   file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
#   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#   save(res, file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
# 
# }
# 
# 
# 
