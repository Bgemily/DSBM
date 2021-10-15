# Description: Estimate the number of clusters using ICL.
# 

apply_ppsbm_ICL = function(### Parameters for generative model
  SEED, N_subj=1, N_node_vec = rep(90,N_subj),
  N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
  total_time=200, 
  conn_patt_var=1, conn_patt_sep = 1.5, const=40, conn_prob_mean = 1, conn_prob_rad = 0, 
  time_shift_struc=max, time_shift_mean_vec = rep(20,N_clus), 
  time_shift_rad = min(time_shift_mean_vec),
  t_vec = seq(0,total_time,length.out=1000),
  ### Parameters for algorithms
  Qmin=N_clus-2, Qmax=N_clus+2 )
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
  
  res = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), d_part=5, 
                Qmin=Qmin, Qmax=Qmax, directed=FALSE, sparse=FALSE, method="hist")
  
  
  # Select best cluster number using ICL ------------------------------------
  
  # ICL-model selection
  data = list(Nijk=Nijk, Time=total_time)
  sol.selec_Q <- modelSelection_Q(data=data,
                                  n=nrow(edge_time_mat),
                                  Qmin=Qmin,
                                  Qmax=Qmax,
                                  directed=FALSE,
                                  sparse=FALSE,
                                  sol.hist.sauv=res)
  
  # best number Q of clusters:
  N_clus_est = sol.selec_Q$Qbest
  
  
  
  
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
              N_clus_est=N_clus_est,
              correct_N_clus=I(N_clus_est==N_clus)*1,
              ICL_vec = sol.selec_Q$all.ICL,
              J_vec = sol.selec_Q$all.J,
              compl_log_lik_vec = sol.selec_Q$all.compl.log.likelihood,
              penalty_vec = sol.selec_Q$all.pen))
  
}



