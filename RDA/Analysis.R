
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


# Get network information for the L/R side ------------------------------------

data_folder = "../Processed_FunctionalData/"
path_vec = list.files(data_folder, full.names = TRUE)
file_vec = list.files(data_folder, full.names = FALSE)

edge_time_mat_list = vector(mode = "list", length = length(path_vec)*2)


for(m in 1:length(path_vec)){ 
  path = path_vec[m]
  
  ### Read information from data
  edge_time_mat_L = as.matrix(read.csv(paste(path, '/EdgeTime_L.csv', sep='')))
  edge_time_mat_L = edge_time_mat_L[,-1]
  edge_time_mat_R = as.matrix(read.csv(paste(path, '/EdgeTime_R.csv', sep='')))
  edge_time_mat_R = edge_time_mat_R[,-1]

  edge_time_mat_list[c(2*m-1,2*m)] = list(edge_time_mat_L, edge_time_mat_R)
}



# # Apply algorithm (ppsbm) ---------------------------------------------------------
# 
# method = "ppsbm"
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
#   res = list(edge_time_mat = edge_time_mat_tmp,
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
# Apply algorithm (our) ---------------------------------------------------------


method = "PDF+pairwise"

N_clus_min = 1 # Number of clusters
N_clus_max = 5
MaxIter = 10 # Maximal iteration number
# bw = 5 # Smoothing bandwidth
freq_trun_vec = c(3,5,7,9) # Cut-off frequency
conv_thres=1e-2
max_iter=10
step_size = 0.5

max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
total_time = max_time + 10
t_vec = seq(0, total_time, 1)

for (subj in 1:length(edge_time_mat_list)) {
  edge_time_mat_list_tmp = edge_time_mat_list[subj]
  ### Get estimation for candidate N_clus and freq_trun
  res_list = list()
  for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
    res_list[[ind_N_clus]] = list()
    N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
    for (ind_freq_trun in 1:length(freq_trun_vec)) {
      freq_trun = freq_trun_vec[ind_freq_trun]
      
      ### Get initialization -----------
      res = get_init_v4(edge_time_mat_list = edge_time_mat_list_tmp, 
                        N_clus = N_clus_tmp, 
                        t_vec = t_vec)
      
      clusters_list_init = res$clusters_list
      n0_vec_list_init = res$n0_vec_list
      n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)

      # Apply algorithm ---------
      
      clusters_list_init -> clusters_list_est
      n0_vec_list_init -> n0_vec_list_est
      n0_mat_list_init -> n0_mat_list_est
      
      time_start = Sys.time()
      ### Estimation z,v,f based on pdf
      res = do_cluster_v14.2.1(edge_time_mat_list = edge_time_mat_list_tmp, 
                               N_clus = N_clus_tmp,
                               clusters_list_init = clusters_list_est,
                               n0_vec_list_init = n0_vec_list_est, 
                               n0_mat_list_init = n0_mat_list_est,
                               total_time = total_time, 
                               max_iter=max_iter, t_vec=t_vec,
                               freq_trun=freq_trun, 
                               step_size=step_size,
                               conv_thres=conv_thres, 
                               MaxIter=MaxIter)
      time_end = Sys.time()
      time_estimation = time_end - time_start
     
      N_iteration = res$N_iteration
      clusters_list_est = res$clusters_list
      v_vec_list_est = res$v_vec_list
      
      ### Get estimated pdf using kernel smoothing
      v_mat_list_est = n0_vec2mat(n0_vec = v_vec_list_est)
      n0_mat_list_est = lapply(v_mat_list_est, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
      center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list_tmp,
                                                     clusters_list = clusters_list_est,
                                                     n0_mat_list = n0_mat_list_est,
                                                     t_vec = t_vec)
      res$center_pdf_array = center_pdf_array_est
      
      # Save results of N_clus_tmp ----------------------------------------------
      res_list[[ind_N_clus]][[ind_freq_trun]] = res
      
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

  ### Retrieve estimation results of the best cluster number
  res = sel_mod_res$res_best
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$center_pdf_array -> center_pdf_array_est
  res$N_iteration -> N_iteration

  ### Save result
  res = list(edge_time_mat = edge_time_mat_list_tmp[[1]],
             clusters_list = clusters_list_est[[1]],
             center_pdf_array = center_pdf_array_est,
             v_vec = v_vec_list_est[[1]],
             N_iteration = N_iteration,
             t_vec = t_vec)
  folder_path = paste0('../Results/Rdata/RDA/', method, '/', file_vec[(subj+1)%/%2])
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  file_name = ifelse(subj%%2==1, yes = "Left", no = "Right")
  now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
  save(res, file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))

}






