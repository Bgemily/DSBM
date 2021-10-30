# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Load simulation result and get network parameters -----------------------

load("../Results/Rdata/SNR_Vis0/main_v5_v5_adap_freq/pr=0.4,n=30,beta=1.3/n/90/N_trial10_20211028_201242.Rdata")
network_param = results[[1]]$network_param

# Generate networks -------------------------------------------------------

network_list = do.call(what = generate_network2_v3, args = network_param)

edge_time_mat_list = network_list$edge_time_mat_list
pdf_true_array = network_list$pdf_true_array
membership_true_list = network_list$membership_true_list
clus_true_list = network_list$clus_true_list
v_true_list = network_list$time_shift_list
order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))


# Set algorithm parameters ------------------------------------------------
N_clus = network_param$N_clus
N_clus_min=N_clus-2 
N_clus_max=N_clus+2

freq_trun=15 
bw=5 
conv_thres=1e-2 
MaxIter=5
jitter_time_rad = 10 
max_iter=10
total_time = network_param$total_time
opt_radius=total_time/2

t_vec = network_param$t_vec
N_subj = network_param$N_subj

# Apply our method --------------------------------------------------------

res_list = list()
for (N_clus_tmp in N_clus_min:N_clus_max) {
  ### Get initialization -----------
  res = get_init_v3(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp, 
                    v_true_list = v_true_list, 
                    jitter_time_rad = 0,
                    t_vec = t_vec)
  
  clusters_list_init = res$clusters_list
  n0_vec_list_init = res$n0_vec_list
  n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
  
  # Apply algorithm ---------
  
  
  ### Estimation z,v,f based on cdf
  res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                        total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                        clusters_list_init = clusters_list_init,
                        n0_vec_list_init = n0_vec_list_init, n0_mat_list_init = n0_mat_list_init,
                        )
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$n0_vec_list -> n0_vec_list_est
  res$n0_mat_list -> n0_mat_list_est
  res$center_cdf_array -> center_cdf_array_est
  
  
  
  ### Estimation z,v,f based on pdf
  res = do_cluster_v14.2.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                           clusters_list_init = clusters_list_est,
                           n0_vec_list_init = n0_vec_list_est, n0_mat_list_init = n0_mat_list_est,
                           total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                           freq_trun=freq_trun, 
                           conv_thres=conv_thres, MaxIter=MaxIter,
                           opt_radius=opt_radius,
                           )
  
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$loss_history -> loss_history
  res$align_time -> align_time
  res$cluster_time -> cluster_time
  res$center_pdf_array -> center_pdf_array_est
  
  # Save results of N_clus_tmp ----------------------------------------------
  
  res_list = c(res_list, list(res))
  
}

save(res_list, network_param,
     file = '../Results/Rdata/SNR_Vis0/main_v5_timeshift_given/pr=0.4,n=90,beta=1.3,one_instance.Rdata')


# Check clustering result -------------------------------------------------

for (tmp in 1:5) {
  print(paste0("When the number of clusters = ", tmp," , the estimated clusters are: "))
  print(res_list[[tmp]]$clusters_list[[1]])
}


# Let time shifts to be zero and re-estimate intensities ------------------

res_list_2 = res_list
for (tmp in 1:5) {
  res_list_2[[tmp]]$v_vec_list[[1]] = res_list_2[[tmp]]$v_vec_list[[1]]*0
  edge_time_mat = edge_time_mat_list[[1]]
  n0_mat_list = list(matrix(0,nrow(edge_time_mat), ncol(edge_time_mat)))
  clusters_list = res_list_2[[tmp]]$clusters_list
  center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                          clusters_list = clusters_list, 
                                          n0_mat_list = n0_mat_list, 
                                          freq_trun = freq_trun,  t_vec = t_vec)
  ### Convert fourier series back to the (smoothed) pdf
  center_pdf_array = res_list_2[[tmp]]$center_pdf_array
  for (q in 1:tmp) {
    for (k in 1:tmp) {
      fft_truncated = center_fft_array[q,k,]
      func = Re(fft(c(tail(fft_truncated, freq_trun+1), 
                      rep(0, length(t_vec)-2*freq_trun-1),
                      head(fft_truncated, freq_trun)), inverse = TRUE))
      center_pdf_array[q,k,] = func
      res_list_2[[tmp]]$freq_trun_mat[q,k] = (sum(fft_truncated!=0)-1)/2
    }
  }
  
  res_list_2[[tmp]]$center_pdf_array = center_pdf_array
}

# Select best cluster number using ICL ------------------------------------

sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list, 
                           N_node_vec = network_param$N_node_vec, 
                           N_clus_min = N_clus_min, 
                           N_clus_max = N_clus_max, 
                           result_list = res_list, 
                           total_time = total_time)

N_clus_est = sel_mod_res$N_clus_est
ICL_vec = sel_mod_res$ICL_vec 
compl_log_lik_vec = sel_mod_res$compl_log_lik_vec 
penalty_vec = sel_mod_res$penalty_vec


sel_mod_res_2 = select_model(edge_time_mat_list = edge_time_mat_list, 
                           N_node_vec = network_param$N_node_vec, 
                           N_clus_min = N_clus_min, 
                           N_clus_max = N_clus_max, 
                           result_list = res_list_2, 
                           total_time = total_time)

N_clus_est_2 = sel_mod_res_2$N_clus_est
ICL_vec_2 = sel_mod_res_2$ICL_vec 
compl_log_lik_vec_2 = sel_mod_res_2$compl_log_lik_vec 
penalty_vec_2 = sel_mod_res_2$penalty_vec


plot(compl_log_lik_vec, type='b', 
     xlab='Number of clusters', 
     ylab = "Log likelihood", 
     main = "Black: v_hat != 0. Red: v_hat == 0",
     ylim=c(-10150,-9980))
lines(compl_log_lik_vec_2, type='b', col=2)

which.max(compl_log_lik_vec_2)

plot(ICL_vec, type='b', 
     xlab='Number of clusters', 
     ylab = "ICL",
     main = "Black: v_hat != 0. Red: v_hat == 0",
     ylim=c(-10300,-10000))
lines(ICL_vec_2, type='b', col=2)


plot(penalty_vec, type='b', 
     xlab='Number of clusters', 
     ylab = "Penalty",
     main = "Black: v_hat != 0. Red: v_hat == 0",
     ylim=c())
lines(penalty_vec_2, type='b', col=2)
