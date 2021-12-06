

### Obtain connecting pattern (pdf of shifted events) for each pair of clusters using multiple subjects
get_center_pdf_array_v2 = function(edge_time_mat_list, clusters_list, 
                                   n0_mat_list=NULL, t_vec=seq(0, 50, 0.05),
                                   ...){  
  time_unit = t_vec[2]-t_vec[1]
  N_subj = length(edge_time_mat_list)
  
  if (length(unique(sapply(clusters_list,length)))>1) {
    stop("The number of clusters is not the same across subjects.")
  }
  N_clus = length(clusters_list[[1]])
  
  pdf_array = array(dim=c(N_clus,N_clus,length(t_vec)))
  
  
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      edge_time_submat_list = list()
      for (m in 1:N_subj) {
        edge_time_submat = edge_time_mat_list[[m]][clusters_list[[m]][[q]], clusters_list[[m]][[l]], drop=F]
        diag(edge_time_submat) = Inf
        
        tau_submat = time_unit * n0_mat_list[[m]][clusters_list[[m]][[q]], clusters_list[[m]][[l]], drop=F]
        edge_time_submat = edge_time_submat - tau_submat
        
        edge_time_submat_list[[m]] = edge_time_submat
      }
      pdf_array[q,l,] = get_pdf_vec(edge_time_vec=unlist(edge_time_submat_list), t_vec=t_vec, ...)
      if(length(unlist(edge_time_submat_list))==0){
        pdf_array[q,l,] = rep(0,length(t_vec))
      }
    }
  }
  
  return(pdf_array)
}

