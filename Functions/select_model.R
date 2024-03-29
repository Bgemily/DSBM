
select_model = function(edge_time_mat_list, N_node_vec, 
                        N_clus_min, N_clus_max, 
                        result_list,
                        total_time)
{
  # Select best cluster number using ICL ------------------------------------
  ICL_mat = compl_log_lik_mat = log_lik_mat = 
    penalty_mat = penalty_2_mat = 
    matrix(nrow = length(result_list[[1]]),
           ncol = length(N_clus_min:N_clus_max))

  for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
    N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]

    for (ind_freq_trun in 1:length(result_list[[ind_N_clus]])) {
      res_tmp = result_list[[ind_N_clus]][[ind_freq_trun]]
      ### Retrieve estimates under current cluster number and current basis number
      clusters_list_tmp = res_tmp$clusters_list
      clus_size_vec = sapply(clusters_list_tmp[[1]], length)
      v_vec_list_tmp = res_tmp$v_vec_list
      v_mat_list_tmp = n0_vec2mat(n0_vec = v_vec_list_tmp)
      center_pdf_array_tmp = res_tmp$center_pdf_array
      freq_trun_mat = res_tmp$freq_trun_mat
      N_basis_mat = res_tmp$N_basis_mat
      if (is.null(N_basis_mat) & !is.null(freq_trun_mat)){
        N_basis_mat = freq_trun_mat
      } else if (is.null(N_basis_mat) & is.null(freq_trun_mat)){
        N_basis_mat = matrix(2,nrow=N_clus_tmp, ncol=N_clus_tmp)
      }
      pi_vec = res_tmp$pi_vec
      if (is.null(pi_vec)){
        pi_vec = clus_size_vec / sum(clus_size_vec)
      }
      tau_mat = res_tmp$tau_mat
      if (is.null(tau_mat)) {
        tau_mat = matrix(0, nrow = N_node_vec[1], ncol = N_clus_tmp)
        for (q in 1:N_clus_tmp) {
          tau_mat[clusters_list_tmp[[1]][[q]], q] = 1
        } 
      }
      
      t_vec = seq(0,total_time,length.out = dim(center_pdf_array_tmp)[3])
      ### Compute log likelihood
      # First term of log likelihood: \sum_{i<j,t_{i,j}<Inf} log( 1-\sum_{q,k} F_{q,k}(T)*tau_{i,q}*tau_{j,k} )
      F_qk_T = apply(center_pdf_array_tmp*(t_vec[2]-t_vec[1]), c(1,2), function(vec)sum(vec[-length(vec)]))
      F_qk_T[is.na(F_qk_T)] = 0
      tau_F_tauT = tau_mat %*% F_qk_T %*% t(tau_mat)
      tmp = edge_time_mat_list[[1]]
      diag(tmp) = NA
      ind_non_edge = which(tmp==Inf)
      log_lik_tmp = (1/2)*sum(log(1-tau_F_tauT)[ind_non_edge])
      
      
      # Second term of log likelihood: \sum_{q,k}{ \sum_{i<j} \log f_{q,k}(t_{i,j}-\max(v_i,v_j))*tau_{i,q}*tau_{j,k} }
      for (q in 1:N_clus_tmp) {
        for (k in 1:N_clus_tmp) {
          log_lik_qk_vec = log(center_pdf_array_tmp[q,k,])
          log_lik_qk_vec[is.na(log_lik_qk_vec)] = 0
          adjst_edge_time_qk_vec = c(unlist(mapply(function(clus_list, mat1, mat2) (mat1-mat2)[clus_list[[q]],clus_list[[k]]], 
                                                   clusters_list_tmp, edge_time_mat_list, v_mat_list_tmp)))
          if(length(adjst_edge_time_qk_vec)==0){
            next
          }
          if (min(adjst_edge_time_qk_vec)<0) {
            adjst_edge_time_qk_vec = adjst_edge_time_qk_vec - min(adjst_edge_time_qk_vec)
          }
          ### Remove adjusted edge times after total_time
          adjst_edge_time_qk_vec[adjst_edge_time_qk_vec>total_time] = Inf
          ### counts: number of edges whose (adjusted) edge time is close to each breakpoint in t_vec
          counts = hist(adjst_edge_time_qk_vec, breaks=t_vec, plot=FALSE, right=FALSE)$counts
          counts = c(counts,0)
          counts = counts / 2 ### Every finite edge time will be counted twice
          
          ### Add compl_log_lik_tmp by \sum_{i,j:z_i=q,z_j=k}{\log f_{q,k}(t_{i,j}-\max(v_i,v_j))}
          ind_tmp = which(counts > 0 & log_lik_qk_vec>-Inf)
          log_lik_tmp = log_lik_tmp + sum(log_lik_qk_vec[ind_tmp]*counts[ind_tmp])
        }
      }
      
      ### Second penalty: \sum_{i}\sum_{q} \tau^{i,q} * \log(\pi_q)
      log_pi_vec = log(pi_vec)
      log_pi_vec[log_pi_vec==-Inf] = 0
      penalty_tmp_2 = -sum(tau_mat %*% log_pi_vec)
      
      ### Compute penalty
      penalty_tmp = (N_clus_tmp-1)/2*sum(log(N_node_vec)) + 
        sum(N_basis_mat[upper.tri(N_basis_mat)],diag(N_basis_mat))/2*sum(log(N_node_vec*(N_node_vec-1)/2))
      
      ### Compute ICL
      ICL_tmp = log_lik_tmp - penalty_tmp - penalty_tmp_2
      compl_log_lik_tmp = log_lik_tmp - penalty_tmp_2
      
      ### Store ICL, loglik, penalties under current cluster number and current basis number
      ICL_mat[ind_freq_trun, ind_N_clus] = ICL_tmp
      compl_log_lik_mat[ind_freq_trun, ind_N_clus] = compl_log_lik_tmp
      log_lik_mat[ind_freq_trun, ind_N_clus] = log_lik_tmp
      penalty_2_mat[ind_freq_trun, ind_N_clus] = penalty_tmp_2
      penalty_mat[ind_freq_trun, ind_N_clus] = penalty_tmp
    }
    
  }
  
  ### Retrieve index of best freq_trun and best cluster number
  ind_best_freq_trun = which(ICL_mat==max(ICL_mat), arr.ind = TRUE)[1,'row']
  ind_best_N_clus = which(ICL_mat==max(ICL_mat), arr.ind = TRUE)[1,'col']
    
  ### Retrieve best ICL,log lik,penalties under best freq_trun
  ICL_vec = ICL_mat[ind_best_freq_trun, ]
  compl_log_lik_vec = compl_log_lik_mat[ind_best_freq_trun, ]
  log_lik_vec = log_lik_mat[ind_best_freq_trun, ]
  penalty_2_vec = penalty_2_mat[ind_best_freq_trun, ]
  penalty_vec = penalty_mat[ind_best_freq_trun, ]
  
  
  N_clus_est = c(N_clus_min:N_clus_max)[ind_best_N_clus]
  
  ### Retrieve estimation results of the best cluster number
  res = result_list[[ind_best_N_clus]][[ind_best_freq_trun]]
  
  # Output ------------------------------------------------------------------
  
  return(list(N_clus_est = N_clus_est, 
              ind_best_N_clus = ind_best_N_clus,
              ind_best_freq_trun = ind_best_freq_trun,
              res_best = res, 
              ICL_mat = ICL_mat,
              compl_log_lik_mat = compl_log_lik_mat, 
              log_lik_mat = log_lik_mat, 
              penalty_2_mat = penalty_2_mat,
              penalty_mat = penalty_mat,
              ICL_vec = ICL_vec, 
              compl_log_lik_vec = compl_log_lik_vec, 
              log_lik_vec = log_lik_vec, 
              penalty_2_vec = penalty_2_vec,
              penalty_vec = penalty_vec,
              counts = counts))
}