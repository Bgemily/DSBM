gener_edge_time_mat = function(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres, t_max=50) { # need to add seed?
  if (!isSymmetric(tau_mat))
    stop('tau_mat should be symmetric.')
  
  N_clus = length(pdfNrdsamp_fun_list)
  N_node = nrow(tau_mat)
  edge_time_mat = matrix(Inf, nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      samples = pdfNrdsamp_fun_list[[q]][[l]]$random(clus_size_vec[q]*clus_size_vec[l]) 
      samples = matrix(samples, clus_size_vec[q], clus_size_vec[l]) 
      edge_time_mat[clus_true[[q]], clus_true[[l]]] = samples
    }
  }
  edge_time_mat[lower.tri(edge_time_mat)] = t(edge_time_mat)[lower.tri(edge_time_mat)] # make it symmetric
  edge_time_mat = edge_time_mat + tau_mat # add time shifts 
  edge_time_mat[pairwise_dist>dist_thres] = Inf # nodes too far away from each other cannot be connected
  
  edge_time_mat[edge_time_mat>t_max] = Inf
  
  if (!isSymmetric(edge_time_mat))
    stop("Constructed edge_time_mat is not symmetric.")
  
  return(edge_time_mat)
}
