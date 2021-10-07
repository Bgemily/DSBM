

# combine subjects and obtain connecting pattern (pdf of shifted events) for each pair of clusters
get_center_pdf_array_mulsubj = function(edge_time_mat_list, clusters_list, n0_vec_list, t_vec=seq(0, 50, 0.05), bw=1){
  time_unit = t_vec[2]-t_vec[1]
  
  N_clus = length(clusters_list[[1]]) 
  
  pdf_array = array(dim=c(N_clus, N_clus, length(t_vec)))
  for (i in 1:N_clus) {
    for (j in 1:N_clus) {
      tmp = function(edge_time_mat, clusters) edge_time_mat[clusters[[i]], clusters[[j]], drop=F]
      edge_time_submat_list = mapply(tmp, edge_time_mat_list, clusters_list, SIMPLIFY = FALSE)
      
      tau_i_vec_list = mapply(function(n0_vec, clusters) time_unit * n0_vec[clusters[[i]]], n0_vec_list, clusters_list, SIMPLIFY = FALSE)
      tau_j_vec_list = mapply(function(n0_vec, clusters) time_unit * n0_vec[clusters[[j]]], n0_vec_list, clusters_list, SIMPLIFY = FALSE)
      
      edge_time_mulsubj_vec_1 = as.vector(unlist(mapply(sweep, edge_time_submat_list, STATS=tau_i_vec_list, MARGIN=1))) # shift according to cluster i
      edge_time_mulsubj_vec_2 = as.vector(unlist(mapply(sweep, edge_time_submat_list, STATS=tau_j_vec_list, MARGIN=2))) # shift according to cluster j   
      
      var_1 = var(edge_time_mulsubj_vec_1[is.finite(edge_time_mulsubj_vec_1)])
      var_2 = var(edge_time_mulsubj_vec_2[is.finite(edge_time_mulsubj_vec_2)])
      
      
      if(var_1<var_2||is.na(var_1)) edge_time_mulsubj_vec = edge_time_mulsubj_vec_1
      else edge_time_mulsubj_vec = edge_time_mulsubj_vec_2
      
      
      pdf_array[i,j,] = get_pdf_vec(edge_time_vec = edge_time_mulsubj_vec, t_vec, bw=bw)
    }
  }
  
  return(pdf_array)
}

