
### Generate network_list, run our algorithm using multiple cluster numbers, and output measurements of errors.
### Based on v3.4.5 
### Add arguments: N_clus_min, N_clus_max
### Add output value: N_clus_est, correct_N_clus, ICL_vec, compl_log_lik_vec, penalty_vec
### Add estimation of cluster number using ICL
### Remove output value: t_vec (this is contained in network_param)


main_v5 = function(### Parameters for generative model
  SEED, N_subj=1, N_node_vec = rep(90,N_subj),
  N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
  total_time=200, 
  conn_patt_var=1, conn_patt_sep = 1.5, const=40, conn_prob_mean = 1, conn_prob_rad = 0, 
  time_shift_struc=max, time_shift_mean_vec = rep(20,N_clus), 
  time_shift_rad = min(time_shift_mean_vec),
  t_vec = seq(0,total_time,length.out=1000),
  ### Parameters for algorithms
  freq_trun=15, bw=5, 
  conv_thres=1e-2, MaxIter=5,
  jitter_time_rad = 10, max_iter=10,
  opt_radius=total_time/2,
  N_clus_min=N_clus-2, N_clus_max=N_clus+2,
  ...)
{
  
  # Generate networks -------------------------------------------------------
  
  network_list = generate_network2_v3(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                                      N_clus = N_clus, clus_size_mat = clus_size_mat,
                                      total_time = total_time, t_vec = t_vec, 
                                      conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
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
  for (N_clus_tmp in N_clus_min:N_clus_max) {
    ### Get initialization -----------
    res = get_init_v4(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp, 
                      t_vec = t_vec)
    
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
    
    ### Estimation z,v,f based on pdf
    res = do_cluster_v14.2.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                             clusters_list_init = clusters_list_est,
                             n0_vec_list_init = n0_vec_list_est, n0_mat_list_init = n0_mat_list_est,
                             total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                             freq_trun=freq_trun, 
                             conv_thres=conv_thres, MaxIter=MaxIter,
                             opt_radius=opt_radius,
                             ...)
    
    
    res$clusters_list -> clusters_list_est
    res$v_vec_list -> v_vec_list_est
    res$loss_history -> loss_history
    res$align_time -> align_time
    res$cluster_time -> cluster_time
    res$center_pdf_array -> center_pdf_array_est
    
    
    # ### Get estimated pdf using kernel smoothing
    # v_mat_list_est = n0_vec2mat(n0_vec = v_vec_list_est)
    # n0_mat_list_est = lapply(v_mat_list_est, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
    # center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
    #                                                clusters_list = clusters_list_est, 
    #                                                n0_mat_list = n0_mat_list_est, 
    #                                                t_vec = t_vec, bw=bw)
    # res$center_pdf_array = center_pdf_array_est

    
    # Save results of N_clus_tmp ----------------------------------------------

    res_list = c(res_list, list(res))
    
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
  penalty_vec = sel_mod_res$penalty_vec
  penalty_2_vec = sel_mod_res$penalty_2_vec
    
  # Retrieve estimation results of the best cluster number ------------------

  res = sel_mod_res$res_best
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$center_pdf_array -> center_pdf_array_est
  res$loss_history -> loss_history
  res$align_time -> align_time
  res$cluster_time -> cluster_time
  
  
  # Compute estimation error ------------------------------------------------
  
  if (N_clus_est == N_clus) {
    # Compute errors of clusters, i.e. Z ------------------------------------
    
    ARI_vec = numeric(length = N_subj)
    for (m in 1:N_subj) {
      ARI_tmp = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est[[m]]), 
                            memb_true_vec = membership_true_list[[m]])
      ARI_vec[m] = ARI_tmp
    }
    weights = N_node_vec / sum(N_node_vec)
    ARI_mean = sum(ARI_vec*weights)
    
    
    
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
    
    
    # Compute error of event rates of all counting processes, i.e. Lambda ----------
    
    ### Warning: Be cautious about the following vectorization - the column order are compatible only for symmetric pdf_true_array.
    ### Compute F
    F_true = matrix(data = pdf_true_array, nrow = N_clus^2, ncol = dim(pdf_true_array)[3])
    F_est = matrix(data = center_pdf_array_est, nrow = N_clus^2, ncol = dim(pdf_true_array)[3])
    
    Lambda_mean_sq_err = 0
    for (m in 1:N_subj) {
      N_node = N_node_vec[m]
      
      ### Compute Z, and hence C
      Z_true = Z_est = matrix(data=0, nrow=N_node, ncol=N_clus)
      for (q in 1:N_clus) {
        Z_true[clus_true_list[[m]][[q]],q] = 1
        Z_est[clusters_list_est[[m]][[q]],q] = 1
      }
      C_true = kronecker(Z_true, Z_true)
      C_est = kronecker(Z_est, Z_est)
      ### id_rm: row index of C_true corresponding to (i,i)'s
      id_rm = which( kronecker(X=seq(N_node),Y=rep(1,N_node)) == 
                       kronecker(Y=seq(N_node),X=rep(1,N_node)) )  
      C_true = C_true[-id_rm, ] 
      C_est = C_est[-id_rm, ] 
      
      ### Compute v in S^v
      v_mat_true = n0_vec2mat(n0_vec = v_true_list[[m]])
      v_mat_vec_true = c(v_mat_true)
      v_mat_vec_true = v_mat_vec_true[-id_rm]
      
      v_mat_est = n0_vec2mat(n0_vec = v_vec_list_est[[m]])
      v_mat_vec_est = c(v_mat_est)
      v_mat_vec_est = v_mat_vec_est[-id_rm]
      
      ### Compute Lambda = S^v \circ (C \times F)
      t_unit = total_time/ncol(F_true)
      Lambda_true = C_true %*% F_true
      Lambda_est = C_est %*% F_est
      for (ij in 1:nrow(Lambda_true)) {
        n0_tmp_true = round(v_mat_vec_true[ij] / t_unit)
        n0_tmp_est = round(v_mat_vec_est[ij] / t_unit)
        Lambda_true[ij, ] = shift_v2(f_origin = Lambda_true[ij, ], n0 = n0_tmp_true)
        Lambda_est[ij, ] = shift_v2(f_origin = Lambda_est[ij, ], n0 = n0_tmp_est)
      }
      if (nrow(Lambda_true)!=(N_node^2-N_node)) {
        stop("The size of Lambda_true is incorrect.")
      }
      
      ### Compute error of Lambda
      error_tmp = sum((Lambda_true-Lambda_est)^2)
      Lambda_mean_sq_err = Lambda_mean_sq_err + error_tmp
    }
    Lambda_mean_sq_err = Lambda_mean_sq_err / sum(N_node_vec^2-N_node_vec)
  }
  else{
    Lambda_mean_sq_err = NA
    ARI_vec=NA 
    ARI_mean=NA
    F_mean_sq_err=NA 
    v_mean_sq_err=NA
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
              # model selection result
              N_clus_est=N_clus_est, 
              correct_N_clus=I(N_clus_est==N_clus)*1, 
              ICL_vec=ICL_vec, 
              compl_log_lik_vec=compl_log_lik_vec, 
              penalty_vec=penalty_vec,
              penalty_2_vec=penalty_2_vec,
              # parameter estimates of best cluster number
              clusters_list_est=clusters_list_est,
              v_vec_list_est=v_vec_list_est,
              # estimation error
              Lambda_mean_sq_err = Lambda_mean_sq_err,
              ARI_vec=ARI_vec, ARI_mean=ARI_mean, 
              F_mean_sq_err=F_mean_sq_err, 
              v_mean_sq_err=v_mean_sq_err,
              # other
              loss_history=loss_history, 
              cluster_time=cluster_time, 
              align_time=align_time
              ))
  
}


