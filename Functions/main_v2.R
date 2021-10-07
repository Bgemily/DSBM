
### Generate network_list, run our algorithm, and output measurements of errors.

main_v2 = function(SEED, N_subj=1, N_node_vec = rep(90,N_subj),
                   N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
                   total_time=200, 
                   conn_patt_var=1, conn_patt_sep = 1.3, const=40, conn_prob_mean = 1, conn_prob_rad = 0, 
                   time_shift_struc=max, time_shift_mean_vec = rep(20,N_clus), time_shift_rad = min(time_shift_mean_vec),
                   max_iter=1,
                   ...)
{
  
# Generate networks -------------------------------------------------------
  
  
  network_list = generate_network2_v2(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                             N_clus = N_clus, clus_size_mat = clus_size_mat,
                             total_time = total_time, 
                             conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                             conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                             time_shift_struc = time_shift_struc,
                             time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)

  edge_time_mat_list = network_list$edge_time_mat_list
  cdf_true_array = network_list$cdf_true_array
  pdf_true_array = network_list$pdf_true_array
  membership_true_list = network_list$membership_true_list
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$time_shift_list
  order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))

  

# Apply algorithm ---------------------------------------------------------

  ### V1
  # res = do_cluster_v2(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus,
  #                     total_time = total_time, max_iter=max_iter, ...)
  
  ### V2
  # print("####################")
  # print("Forcing the initial time shifts to comply with the provided order...")
  # print("####################")
  # res = do_cluster_v2(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus,
  #                     total_time = total_time, max_iter=max_iter,
  #                     order_list = order_true_list, ...)
  
  
  ### V3
  # print("####################")
  # print("Will force the time shifts to comply with the provided order after convergence...")
  # print("####################")
  # res = do_cluster_v3(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus,
  #                     total_time = total_time, max_iter=max_iter,
  #                     order_list = order_true_list, ...)
  
  
  ### V4
  print("####################")
  print("First apply algorithm to each subject. Then combine subjects together.")
  print("####################")
  res = do_cluster_v4(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus,
                      total_time = total_time, max_iter=max_iter,
                      order_list = NULL, conv_thres = 1e-3, ...)
  
  
  ### V5
  # print("####################")
  # print("First apply algorithm to each subject. Then combine subjects together.")
  # print("Normalize node_cdf_array when updating clusters and time shifts.")
  # print("####################")
  # res = do_cluster_v5(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus,
  #                     total_time = total_time, max_iter=max_iter,
  #                     order_list = NULL, ...)
  
  
  
  ### V6
  # print("####################")
  # print("First apply algorithm to each subject. Then combine subjects together.")
  # print("Use pairwise distances when updating clusters.")
  # print("####################")
  # res = do_cluster_v6(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus,
  #                     total_time = total_time, max_iter=max_iter,
  #                     order_list = NULL, ...)
  
  
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$center_cdf_array -> center_cdf_array_est
  
  
# Compute errors of clusters, i.e. Z ------------------------------------
  
  ARI_vec = numeric(length = N_subj)
  for (m in 1:N_subj) {
    ARI_tmp = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est[[m]]), 
                          memb_true_vec = membership_true_list[[m]])
    ARI_vec[m] = ARI_tmp
  }
  ARI_mean = mean(ARI_vec)
  
  
# Compute errors of conn patts, i.e. F ------------------------------------

  ### Match clusters
  res = find_permn(center_cdf_array_from = center_cdf_array_est, 
             center_cdf_array_to = cdf_true_array)
  permn = res$permn
  center_cdf_array_est_permn = center_cdf_array_est[permn, permn, ]
  
  ### Calculate distance
  normed_dist_mat = matrix(nrow=N_clus, ncol=N_clus)
  unnormed_dist_mat = matrix(nrow=N_clus, ncol=N_clus)
  conn_prob_err_mat = matrix(nrow=N_clus, ncol=N_clus)
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      conn_prob_est = tail(center_cdf_array_est_permn[q,k,], 1)
      conn_prob_true = tail(cdf_true_array[q,k,], 1)
      
      cdf_est_normed = center_cdf_array_est_permn[q,k,] / conn_prob_est
      cdf_true_normed = cdf_true_array[q,k,] / conn_prob_true
      
      normed_dist_mat[q,k] = sqrt(sum( (cdf_est_normed - cdf_true_normed)^2 ))
      unnormed_dist_mat[q,k] = sqrt(sum( (center_cdf_array_est_permn[q,k,] - cdf_true_array[q,k,])^2 ))
      conn_prob_err_mat[q,k] = abs(conn_prob_est - conn_prob_true)
    }
  }
  F_mean_sq_err = mean(unnormed_dist_mat[upper.tri(unnormed_dist_mat, diag=TRUE)]^2)
  F_shape_mean_sq_err = mean(normed_dist_mat[upper.tri(normed_dist_mat, diag=TRUE)]^2)
  
  
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
    
  
  
# Compute error of event rates of all counting processes, i.e. Lambda ----------
  
  ### Warning: Be cautious about the following vectorization - the column order are compatible only for symmetric cdf_true_array.
  ### Compute F
  F_true = matrix(data = cdf_true_array, nrow = N_clus^2, ncol = dim(cdf_true_array)[3])
  F_est = matrix(data = center_cdf_array_est, nrow = N_clus^2, ncol = dim(cdf_true_array)[3])
  
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
  
  
  
# Extract network related parameters -----------------------------------------

  network_param = list(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                       N_clus = N_clus, clus_size_mat = clus_size_mat,
                       total_time = total_time, 
                       conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                       conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                       time_shift_struc = time_shift_struc,
                       time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  

# Output ------------------------------------------------------------------

  return(list(network_param=network_param, 
              Lambda_mean_sq_err = Lambda_mean_sq_err,
              ARI_vec=ARI_vec, ARI_mean=ARI_mean,
              normed_dist_mat=normed_dist_mat, unnormed_dist_mat=unnormed_dist_mat, conn_prob_err_mat=conn_prob_err_mat,
              F_mean_sq_err=F_mean_sq_err, F_shape_mean_sq_err=F_shape_mean_sq_err,
              spearman_corr_vec=spearman_corr_vec, spearman_corr_vec_mean=spearman_corr_vec_mean,
              pearson_corr_vec=pearson_corr_vec, pearson_corr_vec_mean=pearson_corr_vec_mean,
              v_mean_sq_err=v_mean_sq_err))
}



# Test --------------------------------------------------------------------

# res = main_v2(SEED=831, N_subj=3, conn_patt_var = 1, conn_patt_sep = 2)
# res$network_param
# res[-1]
