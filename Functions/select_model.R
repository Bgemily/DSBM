
select_model = function(edge_time_mat_list, N_node_vec, 
                         N_clus_min, N_clus_max, 
                         result_list,
                         t_vec, freq_trun_mat=NULL)
{
  # Select best cluster number using ICL ------------------------------------
  ICL_vec = compl_log_lik_vec = penalty_vec = numeric(length = length(N_clus_min:N_clus_max))
  
  for (i in 1:length(ICL_vec)) {
    N_clus_tmp = c(N_clus_min:N_clus_max)[i]
    
    ### Retrieve estimates of i-th cluster number candidate
    clusters_list_tmp = result_list[[i]]$clusters_list
    v_vec_list_tmp = result_list[[i]]$v_vec_list
    v_mat_list_tmp = n0_vec2mat(n0_vec = v_vec_list_tmp)
    center_pdf_array_tmp = result_list[[i]]$center_pdf_array
    freq_trun_mat = result_list[[i]]$freq_trun_mat
    pi_vec = result_list[[i]]$pi_vec
    
    ### Compute log likelihood
    # First term of log likelihood: \sum_{i,j}( -Lambda_{i,j}(T) )
    t_vec = seq(0,200,length.out=dim(center_pdf_array_tmp)[3])
    F_qk_T = apply(center_pdf_array_tmp, c(1,2), sum)*(t_vec[2]-t_vec[1])
    clus_size_vec = sapply(clusters_list_tmp[[1]], length)
    if (length(clus_size_vec)>1) {
      N_nodepair_qk = clus_size_vec %*% t(clus_size_vec) - diag(clus_size_vec)
    }
    else
      N_nodepair_qk = clus_size_vec^2 - clus_size_vec
    
    compl_log_lik_tmp = -sum(F_qk_T*N_nodepair_qk)
    # Set compl_log_lik_tmp to zero in order to align with ppsbm::modelSelection_Q()
    compl_log_lik_tmp = 0
    
    # Second term of log likelihood: \sum_{i,j}{\log f_{z_i,z_j}(t_{i,j}-\max(v_i,v_j))}
    for (q in 1:N_clus_tmp) {
      for (k in q:N_clus_tmp) {
        log_lik_qk_vec = log(center_pdf_array_tmp[q,k,])
        adjst_edge_time_qk_vec = c(unlist(mapply(function(clus_list, mat1, mat2) (mat1-mat2)[clus_list[[q]],clus_list[[k]]], 
                                          clusters_list_tmp, edge_time_mat_list, v_mat_list_tmp)))
        if(length(adjst_edge_time_qk_vec)==0){
          next
        }
        if (min(adjst_edge_time_qk_vec)<0) {
          adjst_edge_time_qk_vec = adjst_edge_time_qk_vec - min(adjst_edge_time_qk_vec)
        }
        ### counts: number of edges whose (adjusted) edge time is close to each breakpoint in t_vec
        counts = hist(adjst_edge_time_qk_vec, breaks=t_vec, plot=FALSE, right=FALSE)$counts
        counts = c(counts,0)
        counts = counts / 2 ### Every finite edge time is counted twice
        
        
        ### Add compl_log_lik_tmp by \sum_{i,j:z_i=q,z_j=k}{\log f_{q,k}(t_{i,j}-\max(v_i,v_j))}
        ind_tmp = which(counts > 0)
        if (sum(log_lik_qk_vec[ind_tmp]*counts[ind_tmp]) == -Inf) {
          next
        }
        compl_log_lik_tmp = compl_log_lik_tmp + sum(log_lik_qk_vec[ind_tmp]*counts[ind_tmp])
      }
    }
    
    ### Third term of log likelihood: \sum_{q} |z^{-1}(q)| * \log(\pi_q)
    if (is.null(pi_vec)){
      pi_vec = clus_size_vec / sum(clus_size_vec)
    }
    compl_log_lik_tmp = compl_log_lik_tmp + sum(clus_size_vec * log(pi_vec))
    
    compl_log_lik_vec[i] = compl_log_lik_tmp
    
    ### Compute penalty
    penalty_tmp = (N_clus_tmp-1)/2*sum(log(N_node_vec)) + 
      sum(freq_trun_mat)/2*sum(log(N_node_vec*(N_node_vec-1)/2))
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