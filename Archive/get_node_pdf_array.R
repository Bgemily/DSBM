
# obtain connecting pattern (pdf of shifted events) between each node and each cluster
get_node_pdf_array = function(edge_time_mat, clusters, n0_vec, n0_mat=NULL, t_vec=seq(0, 50, 0.05), bw=1){  
  time_unit = t_vec[2]-t_vec[1]
  pdf_array = array(dim=c(nrow(edge_time_mat),length(clusters),length(t_vec)))
  for (i in 1:nrow(edge_time_mat)) {
    for (l in 1:length(clusters)) {
      edge_time_submat = edge_time_mat[i, clusters[[l]], drop=F]
      
      if (!is.null(n0_mat)) {
        tau_submat = time_unit * n0_mat[i, clusters[[l]], drop=F]
        edge_time_submat = edge_time_submat - tau_submat
      }
      else{
        tau_vec = time_unit * n0_vec[clusters[[l]]]
        
        edge_time_submat_1 = edge_time_submat
        var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
        
        edge_time_submat_2 = edge_time_submat-tau_vec
        var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
        
        if (var_1<var_2||is.na(var_1)) edge_time_submat = edge_time_submat_1
        else edge_time_submat = edge_time_submat_2
      }
      # pdf = get_pdf_vec(edge_time_submat, t_vec, bw=bw)
      
      
      pdf_array[i,l, ] = get_pdf_vec(edge_time_submat, t_vec, bw=bw)
    }
  }
  
  return(pdf_array)
}

