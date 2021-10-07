

# obtain connecting pattern (pdf of shifted events) for each pair of clusters
get_center_pdf_array = function(edge_time_mat, clusters, n0_vec, n0_mat=NULL, t_vec=seq(0, 50, 0.05), bw=1, clusters_row=clusters){  
  time_unit = t_vec[2]-t_vec[1]
  pdf_array = array(dim=c(length(clusters_row),length(clusters),length(t_vec)))
  
  diag(edge_time_mat) = Inf
  for (q in 1:length(clusters_row)) {
    for (l in 1:length(clusters)) {
      edge_time_submat = edge_time_mat[clusters_row[[q]], clusters[[l]], drop=F]
      
      if (!is.null(n0_mat)) {
        tau_submat = time_unit * n0_mat[clusters_row[[q]], clusters[[l]], drop=F]
        edge_time_submat = edge_time_submat - tau_submat
      }
      else{
        tau_q_vec = time_unit * n0_vec[clusters_row[[q]]]
        tau_l_vec = time_unit * n0_vec[clusters[[l]]]
        
        edge_time_submat_1 = sweep(edge_time_submat, 1, tau_q_vec) # shift according to cluster q
        var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
        
        edge_time_submat_2 = sweep(edge_time_submat, 2, tau_l_vec) # shift according to cluster l
        var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
        
        
        if(var_1<var_2||is.na(var_1)) edge_time_submat = edge_time_submat_1
        else edge_time_submat = edge_time_submat_2
      }
      
      # NEED JUSTIFICATION
      # deal with the situation that shifted event times are negative (thus will be eliminated in estimated pdf)
      # if (min(edge_time_submat) < min(t_vec)) {
      #   print(paste("shifted edge time is negative:",q,l))
      #   MAX_edge_time_submat = max(edge_time_submat[is.finite(edge_time_submat)])
      #   if (MAX_edge_time_submat-min(edge_time_submat) <= max(t_vec)-min(t_vec)) 
      #     edge_time_submat = edge_time_submat - min(edge_time_submat) + min(t_vec)
      #   else{ 
      #     N_early_events = sum( edge_time_submat < MAX_edge_time_submat-(max(t_vec)-min(t_vec)) )
      #     N_late_events = sum( (edge_time_submat > min(edge_time_submat)+max(t_vec)-min(t_vec)) & is.finite(edge_time_submat))
      #     if (N_early_events > N_late_events){ # keep lower tail
      #       edge_time_submat = edge_time_submat - min(edge_time_submat) + min(t_vec)
      #       edge_time_submat[edge_time_submat>max(t_vec)] = Inf
      #     }
      #     else{
      #       edge_time_submat = edge_time_submat + max(t_vec) - MAX_edge_time_submat
      #       edge_time_submat[edge_time_submat<0] = Inf
      #     }
      #   }
      # }
      # 
      
      pdf_array[q,l,] = get_pdf_vec(edge_time_vec = edge_time_submat, t_vec, bw=bw)
    }
  }
  
  return(pdf_array)
}

