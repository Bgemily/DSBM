
# obtain connecting pattern (cdf of shifted events) between each node and each cluster
get_node_cdf_array_v2 = function(edge_time_mat, clusters, n0_mat=0, t_vec=seq(0, 50, 0.05)){  
  time_unit = t_vec[2]-t_vec[1]
  
  adjs_edge_time_mat = edge_time_mat - time_unit*n0_mat
  
  cdf_array = array(dim=c(nrow(edge_time_mat),length(clusters),length(t_vec)))
  for (i in 1:nrow(edge_time_mat)) {
    for (l in 1:length(clusters)) {
      if(length(clusters[[l]])==1 && clusters[[l]]==i){
        cdf_array[i,l, ] = numeric(length(t_vec))
        next
      }
      clus_tmp = setdiff(clusters[[l]],i)
      adjs_edge_time_submat = adjs_edge_time_mat[i, clus_tmp]
      
      
      ############# V1
      # cdf_array[i,l, ] = ecdf(adjs_edge_time_submat)(t_vec)
      ############# V2
      cdf_array[i,l, ] = tryCatch(ecdf(adjs_edge_time_submat)(t_vec),
                                  error=function(x)return(numeric(length(t_vec))))
      ###################################
      
    }
  }
  
  
  
  return(cdf_array)
}

