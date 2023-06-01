
### The default working directory is the current source file location.


# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "../Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(tidyverse)
library(scales)
library(combinat)
library(cluster)
library(mclust)
library(reshape2)


# Get network information for the L/R side ------------------------------------

data_folder = "../Processed_FunctionalData/"
path_vec = list.files(data_folder, full.names = TRUE)
file_vec = list.files(data_folder, full.names = FALSE)
edge_time_mat_list = vector(mode = "list", length = length(path_vec)*2)
avai_inds_list = list()
for(m in 1:length(path_vec)){ 
  path = path_vec[m]
  
  ### Read information from data
  edge_time_mat = as.matrix(read.csv(paste(path, '/EdgeTime.csv', sep='')))
  edge_time_mat = edge_time_mat[,-1]
  avai.inds = as.matrix(read.csv(paste(path,'/AvaiNeurons.csv',sep='')))
  avai.inds = avai.inds[,-1];
  locs.all = as.matrix(read.csv(paste(path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs_mat = locs.all[avai.inds,]
  ### Separate networks in left and right spines
  inds_L = which(locs_mat[,2]<0)
  inds_R = which(locs_mat[,2]>0)
  locs_mat_L = locs_mat[inds_L,]
  locs_mat_R = locs_mat[inds_R,]
  edge_time_mat_L = edge_time_mat[inds_L, inds_L]
  edge_time_mat_R = edge_time_mat[inds_R, inds_R]
  ### Remove neurons with no connection and neurons with abnormal behavior
  avai_inds_L = which(rowSums(edge_time_mat_L<Inf)>0 )
  avai_inds_R = which(rowSums(edge_time_mat_R<Inf)>0 )
  if (m==2) {
    avai_inds_L = avai_inds_L[-c(1,2)]
    avai_inds_R = avai_inds_R[-c(45,34,15)]
  }
  edge_time_mat_L = edge_time_mat_L[avai_inds_L, avai_inds_L]
  edge_time_mat_R = edge_time_mat_R[avai_inds_R, avai_inds_R]
  
  edge_time_mat_list[c(2*m-1,2*m)] = list(edge_time_mat_L, edge_time_mat_R)
  avai_inds_list[c(2*m-1,2*m)] = list(inds_L[avai_inds_L], inds_R[avai_inds_R])
  
}






# Apply algorithm ---------------------------------------------------------
N_clus_min = 1 # Minimal number of clusters
N_clus_max = 4 # Maximal number of clusters
MaxIter = 10 # Maximal iteration number between centering step and clustering step
conv_thres = 1e-3 # Convergence threshold
max_iter = 10 # Maximal iteration number for within centering step
N_restart = 1 #

max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
total_time = max_time + 10
total_time = 336
t_vec = seq(0, total_time, 1)


for (subj in 3:4){
  ICL_best = -Inf
  res_best = NULL
  res_multi_restart = list()
  ICL_history = matrix(nrow=N_restart,ncol=length(N_clus_min:N_clus_max))
  compl_log_lik_history = matrix(nrow=N_restart,ncol=length(N_clus_min:N_clus_max))
  
  time_start = Sys.time()
  for (ind_restart in 1:N_restart){
    ### Apply method with freq_trun and total_time_cutoff
    edge_time_mat_list_tmp = edge_time_mat_list[subj]
    avai_inds = avai_inds_list[[subj]]
    ### Get estimation for candidate N_clus and freq_trun
    res_list = list()
    for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
      res_list[[ind_N_clus]] = list()
      N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
      for (ind_freq_trun_tmp in 1:1) {
        ### Get initialization
        res = get_init(edge_time_mat_list = edge_time_mat_list_tmp,
                       N_clus = N_clus_tmp,
                       t_vec = t_vec)
        
        clusters_list_init = res$clusters_list
        n0_vec_list_init = res$n0_vec_list
        n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
        
        # Apply algorithm
        ### Estimation z,v,f based on cdf
        time_start = Sys.time()
        res = do_cluster_cdf(edge_time_mat_list = edge_time_mat_list_tmp,
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
        res$t_vec=t_vec
        
        res$edge_time_mat = edge_time_mat_list_tmp[[1]]
        res$avai_inds = avai_inds
        res$v_vec = res$v_vec_list[[1]]
        
        # Save results of N_clus_tmp
        res_list[[ind_N_clus]][[ind_freq_trun_tmp]] = res
        
      }
    }
    
    ### Select best cluster number using ICL
    sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list_tmp,
                               N_node_vec = sapply(edge_time_mat_list_tmp,nrow),
                               N_clus_min = N_clus_min,
                               N_clus_max = N_clus_max,
                               result_list = res_list,
                               total_time = total_time)
    
    N_clus_est = sel_mod_res$N_clus_est
    ICL_vec = sel_mod_res$ICL_vec
    compl_log_lik_vec = sel_mod_res$compl_log_lik_vec
    log_lik_vec = sel_mod_res$log_lik_vec
    penalty_2_vec = sel_mod_res$penalty_2_vec
    penalty_vec = sel_mod_res$penalty_vec
    ICL_mat = sel_mod_res$ICL_mat
    compl_log_lik_mat = sel_mod_res$compl_log_lik_mat
    log_lik_mat = sel_mod_res$log_lik_mat
    penalty_2_mat = sel_mod_res$penalty_2_mat
    penalty_mat = sel_mod_res$penalty_mat
    
    ### Retrieve estimation results of the best cluster number
    res = sel_mod_res$res_best
    res$clusters_list -> clusters_list_est
    res$v_vec_list -> v_vec_list_est
    res$center_pdf_array -> center_pdf_array_est
    res$N_iteration -> N_iteration
    
    ### Save result
    res = list(res_list = res_list,
               edge_time_mat = edge_time_mat_list_tmp[[1]],
               avai_inds = avai_inds,
               clusters_list = clusters_list_est[[1]],
               center_pdf_array = center_pdf_array_est,
               v_vec = v_vec_list_est[[1]],
               N_clus_est = N_clus_est,
               ICL_vec = ICL_vec,
               compl_log_lik_vec = compl_log_lik_vec,
               log_lik_vec = log_lik_vec,
               penalty_2_vec = penalty_2_vec,
               penalty_vec = penalty_vec,
               ICL_mat = ICL_mat,
               compl_log_lik_mat = compl_log_lik_mat,
               log_lik_mat = log_lik_mat,
               penalty_2_mat = penalty_2_mat,
               penalty_mat = penalty_mat,
               N_iteration = N_iteration,
               t_vec = t_vec)
    
    ### Update res_best and ICL_best
    if (max(res$ICL_vec)>ICL_best){
      ICL_best = max(res$ICL_vec)
      res_best = res
    }
    ICL_history[ind_restart, ] = res$ICL_vec
    compl_log_lik_history[ind_restart, ] = res$compl_log_lik_vec
    res_multi_restart[ind_restart] = list(res)
  }
  time_end = Sys.time()
  time_total = time_end - time_start
  
  ### Save result with freq_trun and total_time_cutoff
  res = res_best
  res_Nclus = res$res_list
  method = paste0("CDF_v15_rmv2+3_keeptop",
                  "_Nrestart",N_restart,
                  "_","totaltime",total_time)
  folder_path = paste0('../Results/Rdata/RDA_v3/', method, '/', file_vec[(subj+1)%/%2])
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
  now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
  save(res,ICL_history,compl_log_lik_history,time_total,
       res_Nclus,
       res_multi_restart,
       file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
  
}


