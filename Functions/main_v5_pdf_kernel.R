
### Generate network_list, run our algorithm using multiple cluster numbers, and output measurements of errors.
### Based on v3.4.5 
### Add arguments: N_clus_min, N_clus_max
### Add output value: N_clus_est, correct_N_clus, ICL_vec, compl_log_lik_vec, penalty_vec
### Add estimation of cluster number using ICL
### Remove output value: t_vec (this is contained in network_param)


main_v5_pdf_kernel = function(### Parameters for generative model
  SEED, N_subj=1, N_node_vec = rep(90,N_subj),
  N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
  total_time=200, 
  conn_patt_var=1, conn_patt_sep = 1.5, const=20, conn_prob_mean = 1, conn_prob_rad = 0, 
  time_shift_struc=max, time_shift_mean_vec = rep(20,N_clus), 
  time_shift_rad = min(time_shift_mean_vec),
  t_vec = seq(0,total_time,length.out=1000),
  ### Parameters for algorithms
  freq_trun=15, bw=5, 
  conv_thres=1e-2, MaxIter=5,
  gamma = 0.001,
  jitter_time_rad = 10, max_iter=50,
  opt_radius=total_time/2,
  N_clus_min=N_clus-2, N_clus_max=N_clus+2,
  freq_trun_vec=NULL,
  fix_timeshift=FALSE,
  save_center_pdf_array=FALSE,
  ...)
{
  
  # Generate networks -------------------------------------------------------
  
  network_list = generate_network(SEED = SEED, N_node = N_node_vec[1],
                                  N_clus = N_clus, clus_size_vec = clus_size_mat[1,],
                                  total_time = total_time, t_vec = t_vec, 
                                  conn_patt_sep = conn_patt_sep, const = const,
                                  conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                                  time_shift_struc = time_shift_struc,
                                  time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  
  edge_time_mat_list = network_list$edge_time_mat_list
  pdf_true_array = network_list$pdf_true_array
  membership_true_list = network_list$membership_true_list
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$time_shift_list
  order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))
  
  
  # Fit model for various cluster number ------------------------------------
  
  res_list = list()
  if(is.null(freq_trun_vec)){
    freq_trun_vec = c(freq_trun)
  }
  for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
    res_list[[ind_N_clus]] = list()
    N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
    
    for (ind_freq_trun in 1:length(freq_trun_vec)) {
      freq_trun = freq_trun_vec[ind_freq_trun]
      
      ### Get initialization -----------
      if (fix_timeshift) {
        res = get_init_jittertruev(edge_time_mat_list = edge_time_mat_list, 
                                   N_clus = N_clus_tmp, 
                                   v_true_list = v_true_list, 
                                   jitter_time_rad = 0,
                                   t_vec = t_vec)
      } else{
        res = get_init(edge_time_mat_list = edge_time_mat_list, 
                       N_clus = N_clus_tmp, 
                       t_vec = t_vec)
      }
      
      clusters_list_init = res$clusters_list
      n0_vec_list_init = res$n0_vec_list
      n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
      
      ### Evaluate accuracy of initial time shifts
      spearman_corr_vec = numeric(length = N_subj)
      for (m in 1:N_subj) {
        spearman_corr_vec[m] = cor(v_true_list[[m]], n0_vec_list_init[[m]], method = "spearman")
      }
      spearman_corr_vec_mean_init = mean(spearman_corr_vec)
      
      
      
      # Apply algorithm ---------
      
      clusters_list_init -> clusters_list_est
      n0_vec_list_init -> n0_vec_list_est
      n0_mat_list_init -> n0_mat_list_est
      
      time_start = Sys.time()
      ### Estimation z,v,f based on pdf
      res = do_cluster_pdf(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                           clusters_list_init = clusters_list_est,
                           n0_vec_list_init = n0_vec_list_est, n0_mat_list_init = n0_mat_list_est,
                           total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                           freq_trun=freq_trun, 
                           conv_thres=conv_thres, MaxIter=MaxIter,
                           opt_radius=opt_radius,
                           fix_timeshift=fix_timeshift,
                           gamma = gamma,
                           ...)
      time_end = Sys.time()
      time_estimation = time_end - time_start
      time_estimation = as.numeric(time_estimation, units='secs')
      
      res$clusters_list -> clusters_list_est
      res$v_vec_list -> v_vec_list_est
      res$loss_history -> loss_history
      res$align_time -> align_time
      res$cluster_time -> cluster_time
      
      N_iteration = res$N_iteration
      
      
      ### Get estimated pdf using kernel smoothing
      v_mat_list_est = n0_vec2mat(n0_vec = v_vec_list_est)
      n0_mat_list_est = lapply(v_mat_list_est, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
      center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list,
                                                     clusters_list = clusters_list_est,
                                                     n0_mat_list = n0_mat_list_est,
                                                     t_vec = t_vec)
      res$center_pdf_array_kernel = center_pdf_array_est
      
      
      # Save results of N_clus_tmp ----------------------------------------------
      
      res_list[[ind_N_clus]][[ind_freq_trun]] = res
      
    }
  }
  
  
  
  # Select best cluster number using ICL ------------------------------------
  
  sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list, 
                             N_node_vec = N_node_vec, 
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
  
  # Retrieve estimation results of the best cluster number ------------------
  
  res = sel_mod_res$res_best
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$center_pdf_array_kernel -> center_pdf_array_est
  res$loss_history -> loss_history
  res$align_time -> align_time
  res$cluster_time -> cluster_time
  
  
  # Compute estimation error ------------------------------------------------
  
  if (N_clus_est == N_clus) {
    # Compute errors of conn patts, i.e. F ------------------------------------
    
    ### Match clusters [WARNING: This dooes not work for multiple subjects.]
    res = find_permn(center_cdf_array_from = center_pdf_array_est,
                     center_cdf_array_to = pdf_true_array)
    permn = res$permn
    center_pdf_array_est_permn = center_pdf_array_est[permn, permn, ]
    
    
    ### Calculate distance 
    dist_mat = matrix(nrow=N_clus, ncol=N_clus)
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        dist_mat[q,k] = sqrt(sum( (center_pdf_array_est_permn[q,k,] - pdf_true_array[q,k,])^2 * 
                                    (t_vec[2]-t_vec[1]) ))
      }
    }
    F_mean_sq_err = mean(dist_mat[upper.tri(dist_mat, diag=TRUE)]^2)
    
    
    # Compute errors of clusters, i.e. Z ------------------------------------
    
    ARI_vec = numeric(length = N_subj)
    clusters_list_est_permn = clusters_list_est
    for (m in 1:N_subj) {
      ARI_tmp = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est[[m]]), 
                            memb_true_vec = membership_true_list[[m]])
      ARI_vec[m] = ARI_tmp
      
      clusters_list_est_permn[[m]] = clusters_list_est[[m]][permn]
    }
    weights = N_node_vec / sum(N_node_vec)
    ARI_mean = sum(ARI_vec*weights)
    
    
    
    # Compute error of time shifts, i.e. v --------------------------------------------
    
    spearman_corr_vec = numeric(length = N_subj)
    pearson_corr_vec = numeric(length = N_subj)
    for (m in 1:N_subj) {
      spearman_corr_vec[m] = cor(v_true_list[[m]], v_vec_list_est[[m]], method = "spearman")
      pearson_corr_vec[m] = cor(v_true_list[[m]], v_vec_list_est[[m]], method = "pearson")
    }
    spearman_corr_vec_mean = mean(spearman_corr_vec)
    pearson_corr_vec_mean = mean(pearson_corr_vec)
    
    v_mean_sq_err = mean((unlist(v_true_list)-unlist(v_vec_list_est))^2) 
    
    ### Debug
    # n0_mat_list_est = lapply(v_vec_list_est, function(v_vec)n0_vec2mat(n0_vec = v_vec/t_vec[2]))
    # center_pdf_array_est = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
    #                                                clusters_list = clusters_list_est, 
    #                                                t_vec = seq(0,200,length.out=1000), 
    #                                                n0_mat_list = n0_mat_list_est)
    
    
  }
  else{
    ARI_vec=NA 
    ARI_mean=NA
    F_mean_sq_err=NA 
    v_mean_sq_err=NA
    clusters_list_est_permn=NA
  }
  
  # Extract network related parameters -----------------------------------------
  
  network_param = list(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                       N_clus = N_clus, clus_size_mat = clus_size_mat,
                       total_time = total_time, t_vec = t_vec,
                       conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                       conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                       time_shift_struc = time_shift_struc,
                       time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  freqs = N_node_vec
  freqs = freqs/sum(freqs) 
  N_node_vec_entropy = -sum(ifelse(freqs>0, freqs*log(freqs),0)) / log(N_subj)
  
  
  # Output ------------------------------------------------------------------
  
  return(list(network_param=network_param, 
              N_node_vec_entropy=N_node_vec_entropy,
              freq_trun_vec=freq_trun_vec,
              # model selection result
              N_clus_est=N_clus_est, 
              correct_N_clus=I(N_clus_est==N_clus)*1, 
              ICL_vec=ICL_vec, 
              compl_log_lik_vec=compl_log_lik_vec, 
              log_lik_vec=log_lik_vec,
              penalty_2_vec=penalty_2_vec,
              penalty_vec=penalty_vec,
              ICL_mat = ICL_mat,
              compl_log_lik_mat = compl_log_lik_mat, 
              log_lik_mat = log_lik_mat, 
              penalty_2_mat = penalty_2_mat,
              penalty_mat = penalty_mat,
              # parameter estimates of best cluster number
              clusters_list_est=clusters_list_est,
              v_vec_list_est=v_vec_list_est,
              clusters_list_est_permn=clusters_list_est_permn,
              center_pdf_array_est_permn=switch(save_center_pdf_array, 
                                                "TRUE"=center_pdf_array_est_permn,
                                                "FALSE"=NULL),
              # estimation error
              ARI_vec=ARI_vec, ARI_mean=ARI_mean, 
              F_mean_sq_err=F_mean_sq_err, 
              v_mean_sq_err=v_mean_sq_err,
              # other
              time_estimation=time_estimation,
              N_iteration=N_iteration,
              loss_history=loss_history, 
              cluster_time=cluster_time, 
              align_time=align_time
  ))
  
}


