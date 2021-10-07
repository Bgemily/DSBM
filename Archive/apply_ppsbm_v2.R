
### Use generate_network2_v3
### Output measurements of errors
apply_ppsbm_v2 = function(### Parameters for generative model
                          SEED, N_subj=1, N_node_vec = rep(90,N_subj),
                          N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
                          total_time=200, 
                          conn_patt_var=1, conn_patt_sep = 1.5, const=40, conn_prob_mean = 1, conn_prob_rad = 0, 
                          time_shift_struc=max, time_shift_mean_vec = rep(20,N_clus), 
                          time_shift_rad = min(time_shift_mean_vec),
                          ### Parameters for algorithms
                          Qmin=N_clus, Qmax=N_clus, bw=5,
                          t_vec = seq(0,total_time,length.out=1000) )
{
  
  # Generate networks -------------------------------------------------------
  
  
  network_list = generate_network2_v3(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                                      N_clus = N_clus, clus_size_mat = clus_size_mat,
                                      total_time = total_time, 
                                      conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                                      conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                                      time_shift_struc = time_shift_struc,
                                      t_vec = t_vec,
                                      time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  
  edge_time_mat_list = network_list$edge_time_mat_list
  cdf_true_array = network_list$cdf_true_array
  pdf_true_array = network_list$pdf_true_array
  membership_true_list = network_list$membership_true_list
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$time_shift_list
  order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))
  
  
  edge_time_mat = edge_time_mat_list[[1]]
  
  

  # Apply PPSBM -------------------------------------------------------------

  library(ppsbm)
  time.seq = numeric(sum(edge_time_mat<Inf))
  type.seq = numeric(sum(edge_time_mat<Inf))
  current_ind = 1
  for (i in 1:nrow(edge_time_mat)) {
    for (j in 1:ncol(edge_time_mat)) {
      if (edge_time_mat[i,j]<Inf){
        time.seq[current_ind] = edge_time_mat[i,j]
        type.seq[current_ind] = convertNodePair(i, j, n = nrow(edge_time_mat), directed = FALSE)
        current_ind = current_ind+1
      }
    }
  }
  data = list(time.seq=time.seq, type.seq=type.seq, Time=total_time)
  Nijk = statistics(data, nrow(edge_time_mat), K=2^6, directed = FALSE)
  
  # res = mainVEM(data=data, n=nrow(edge_time_mat), Qmin=3, directed=FALSE, method="kernel",
  #               d_part=0, n_perturb=0)
  res = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), d_part=5, 
                Qmin=Qmin, Qmax=Qmax, directed=FALSE, method="hist")[[1]]
  res$clusters = mem2clus(apply(res$tau, 2, which.max)) 
  
  
  clusters_list_est = list(res$clusters)
  

  center_cdf_array_est = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                 clusters_list = clusters_list_est, 
                                                 n0_mat_list = list(edge_time_mat*0), 
                                                 t_vec = t_vec)
  
  center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                 clusters_list = clusters_list_est, 
                                                 n0_mat_list = list(matrix(0,nrow(edge_time_mat), 
                                                                           ncol(edge_time_mat))), 
                                                 t_vec = t_vec, bw=bw)
  
  
  
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
  membership_est_vec = clus2mem(clusters_list_est[[1]])
  z_hat = t(dummies::dummy(membership_est_vec))
  z_true = t(dummies::dummy(membership_true_list[[1]]))
  permn = ppsbm::permuteZEst(z = z_true, hat.z = z_hat)
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
  
  
  # Extract network related parameters -----------------------------------------
  
  network_param = list(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                       N_clus = N_clus, clus_size_mat = clus_size_mat,
                       total_time = total_time, t_vec = t_vec,
                       conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                       conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                       time_shift_struc = time_shift_struc,
                       time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  
  
  # Output ------------------------------------------------------------------
  
  return(list(network_param=network_param, 
              t_vec=t_vec,
              clusters_list_est=clusters_list_est,
              ARI_vec=ARI_vec, ARI_mean=ARI_mean,
              F_mean_sq_err=F_mean_sq_err))
}




# Test --------------------------------------------------------------------

# res = apply_ppsbm_v2(SEED=1983, conn_prob_mean = 0.7, N_node_vec = c(60),
#                      conn_patt_sep = 1.8, time_shift_mean_vec = rep(20,3),
#                      t_vec = seq(0,200,length.out=1000))
# 
#   network_list = generate_network2_v3(SEED = res$network_param$SEED,
#                                       N_subj = res$network_param$N_subj,
#                                       N_node_vec = res$network_param$N_node_vec,
#                                       N_clus = res$network_param$N_clus,
#                                       clus_size_mat = res$network_param$clus_size_mat,
#                                       total_time = res$network_param$total_time,
#                                       t_vec = res$network_param$t_vec,
#                                       conn_patt_var = res$network_param$conn_patt_var,
#                                       conn_patt_sep = res$network_param$conn_patt_sep,
#                                       const = res$network_param$const,
#                                       conn_prob_mean = res$network_param$conn_prob_mean,
#                                       conn_prob_rad = res$network_param$conn_prob_rad,
#                                       time_shift_struc = res$network_param$time_shift_struc,
#                                       time_shift_mean_vec = res$network_param$time_shift_mean_vec,
#                                       time_shift_rad = res$network_param$time_shift_rad)
#   edge_time_mat_list = network_list$edge_time_mat_list
#   membership_true_list = network_list$membership_true_list
#   clusters_list_est = res$clusters_list_est
#   pdf_true_array = res$pdf_true_array
#   t_vec = res$t_vec
#   center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list,
#                                                  clusters_list = clusters_list_est,
#                                                  n0_mat_list = lapply(edge_time_mat_list, 
#                                                                       function(mat)matrix(0,nrow(mat),ncol(mat))),
#                                                  t_vec = t_vec, bw=5)
#   ### Match clusters
#   membership_est_vec = clus2mem(clusters_list_est[[1]])
#   z_hat = t(dummies::dummy(membership_est_vec))
#   z_true = t(dummies::dummy(membership_true_list[[1]]))
#   permn = ppsbm::permuteZEst(z = z_true, hat.z = z_hat)
#   center_pdf_array_est_permn = center_pdf_array_est[permn, permn, ]
# 
#   plot_pdf_array(pdf_array_list = list(center_pdf_array_est_permn),
#                  pdf_true_array = pdf_true_array, t_vec = t_vec)
# 
#   res$ARI_vec
