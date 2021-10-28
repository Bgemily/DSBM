
select_model = function(edge_time_mat_list, N_node_vec, 
                        N_clus_min, N_clus_max, 
                        result_list,
                        total_time)
{
  # Select best cluster number using ICL ------------------------------------
  ICL_vec = compl_log_lik_vec = penalty_vec = numeric(length = length(N_clus_min:N_clus_max))
  
  for (i in 1:length(ICL_vec)) {
    N_clus_tmp = c(N_clus_min:N_clus_max)[i]
    
    ### Retrieve estimates of i-th cluster number candidate
    clusters_list_tmp = result_list[[i]]$clusters_list
    clus_size_vec = sapply(clusters_list_tmp[[1]], length)
    v_vec_list_tmp = result_list[[i]]$v_vec_list
    v_mat_list_tmp = n0_vec2mat(n0_vec = v_vec_list_tmp)
    center_pdf_array_tmp = result_list[[i]]$center_pdf_array
    freq_trun_mat = result_list[[i]]$freq_trun_mat
    N_basis_mat = result_list[[i]]$N_basis_mat
    if (is.null(N_basis_mat)){
      N_basis_mat = 1+2*freq_trun_mat
    }
    pi_vec = result_list[[i]]$pi_vec
    if (is.null(pi_vec)){
      pi_vec = clus_size_vec / sum(clus_size_vec)
    }
    tau_mat = result_list[[i]]$tau_mat
    if (is.null(tau_mat)) {
      tau_mat = matrix(0, nrow = N_node_vec[1], ncol = N_clus_tmp)
      for (q in 1:N_clus_tmp) {
        tau_mat[clusters_list_tmp[[1]][[q]], q] = 1
      } 
    }
    
    t_vec = seq(0,total_time,length.out = dim(center_pdf_array_tmp)[3])
    ### Compute log likelihood
    # First term of log likelihood: -\sum_{i<j}( \sum_{q,k} F_{q,k}(T)*tau_{i,q}*tau_{j,k} )
    F_qk_T = apply(center_pdf_array_tmp*(t_vec[2]-t_vec[1]), c(1,2), function(vec)sum(vec[-length(vec)]))
    tau_F_tauT = tau_mat %*% F_qk_T %*% t(tau_mat)
    compl_log_lik_tmp = -sum(tau_F_tauT[upper.tri(tau_F_tauT)])
    # Set compl_log_lik_tmp to zero in order to align with ppsbm::modelSelection_Q()
    # compl_log_lik_tmp = 0
    
    
    # Second term of log likelihood: \sum_{q,k}{ \sum_{i<j} \log f_{q,k}(t_{i,j}-\max(v_i,v_j))*tau_{i,q}*tau_{j,k} }
    adjst_edge_time_mat = edge_time_mat_list[[1]] - v_mat_list_tmp[[1]]
    if (min(adjst_edge_time_mat)<0) {
      adjst_edge_time_mat = adjst_edge_time_mat - min(adjst_edge_time_mat)
    }               
    counts_array = array(0, dim=c(N_node_vec[1],N_node_vec[1],length(t_vec)))
    for (ii in 1:N_node_vec[1]){
      for (jj in 1:N_node_vec[1]){
        counts = hist(adjst_edge_time_mat[ii,jj], breaks=t_vec, plot=FALSE, right=FALSE)$counts
        counts_array[ii,jj,-length(t_vec)] = counts
      }
    }
    
    for (q in 1:N_clus_tmp) {
      for (k in 1:N_clus_tmp) {
        ### Calculate $\log f_{q,k}(t_{i,j}-\max(v_i,v_j))*tau_{i,q}*tau_{j,k}$ for all (i,j)
        log_lik_qk_vec = log(center_pdf_array_tmp[q,k,])
        log_lik_qk_vec[log_lik_qk_vec==-Inf] = 0
        log_lik_qk_mat = matrix(0, nrow=N_node_vec[1], ncol=N_node_vec[1])               
        for (ii in 1:N_node_vec[1]){
          for (jj in 1:N_node_vec[1]){
            ind_tmp = which(counts_array[ii,jj,]>0)
            if (length(ind_tmp)==0)
              next
            log_lik_qk_mat[ii,jj] = sum(counts_array[ii,jj,ind_tmp]*log_lik_qk_vec[ind_tmp])
          }
        }
        mat_tmp = log_lik_qk_mat * (tau_mat[,q]%*%t(tau_mat[,k])) 
        
        
        ### Add compl_log_lik_tmp by \sum_{i<j} \log f_{q,k}(t_{i,j}-\max(v_i,v_j))*tau_{i,q}*tau_{j,k}        
        compl_log_lik_tmp = compl_log_lik_tmp + sum(mat_tmp[upper.tri(mat_tmp)])
        
      }
    }
    
    ### Third term of log likelihood: \sum_{i}\sum_{q} \tau^{i,q} * \log(\pi_q)
    compl_log_lik_tmp = compl_log_lik_tmp + sum(tau_mat %*% log(pi_vec))
    
    compl_log_lik_vec[i] = compl_log_lik_tmp
    
    ### Compute penalty
    penalty_tmp = (N_clus_tmp-1)/2*sum(log(N_node_vec)) + 
      sum(N_basis_mat[upper.tri(N_basis_mat)],diag(N_basis_mat))/2*sum(log(N_node_vec*(N_node_vec-1)/2))
    penalty_vec[i] = penalty_tmp
    
    ### Compute ICL
    ICL_vec[i] = compl_log_lik_vec[i] - penalty_vec[i]
    
  }
  
  
  ind_best_N_clus = which.max(ICL_vec)
  N_clus_est = c(N_clus_min:N_clus_max)[ind_best_N_clus]
  
  
  ### Retrieve estimation results of the best cluster number
  res = result_list[[ind_best_N_clus]]
  
  # Output ------------------------------------------------------------------
  
  return(list(N_clus_est = N_clus_est, 
              res_best = res, 
              ICL_vec = ICL_vec, 
              compl_log_lik_vec = compl_log_lik_vec, 
              penalty_vec = penalty_vec,
              counts = counts))
}