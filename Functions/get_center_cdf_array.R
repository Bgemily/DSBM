

# obtain connecting pattern (pdf of shifted events) for each pair of clusters
get_center_cdf_array = function(edge_time_mat, clusters, n0_mat=0, t_vec=seq(0, 50, 0.05), clusters_row=clusters, standardize=FALSE){  
  time_unit = t_vec[2]-t_vec[1]
  
  adjs_edge_time_mat = edge_time_mat - time_unit*n0_mat
  diag(adjs_edge_time_mat) = NA # remove (i,i)
  
  cdf_array = array(dim=c(length(clusters_row),length(clusters),length(t_vec)))
  for (q in 1:length(clusters_row)) {
    for (l in 1:length(clusters)) {
      if(length(clusters[[l]])==1){
        cdf_array[q,l, ] = numeric(length(t_vec))
        next
      }
      
      adjs_edge_time_submat = adjs_edge_time_mat[clusters_row[[q]], clusters[[l]]]
      
      # NEED JUSTIFICATION
      # deal with the situation that shifted event times are negative (thus will be eliminated in estimated cdf)
      if (min(adjs_edge_time_submat, na.rm=T) < min(t_vec)) {
        MIN_edge_time_submat = min(adjs_edge_time_submat, na.rm=T)
        MAX_edge_time_submat = max(adjs_edge_time_submat[is.finite(adjs_edge_time_submat)], na.rm=T)
        if (MAX_edge_time_submat-MIN_edge_time_submat <= max(t_vec)-min(t_vec)) 
          adjs_edge_time_submat = adjs_edge_time_submat - MIN_edge_time_submat + min(t_vec)
        else{ 
          N_early_events = sum( isTRUE(adjs_edge_time_submat < MAX_edge_time_submat-(max(t_vec)-min(t_vec))) )
          N_late_events = sum( (adjs_edge_time_submat > MIN_edge_time_submat+max(t_vec)-min(t_vec)) & is.finite(adjs_edge_time_submat))
          if (N_early_events > N_late_events){ # keep lower tail
            adjs_edge_time_submat = adjs_edge_time_submat - MIN_edge_time_submat + min(t_vec)
            adjs_edge_time_submat[which(adjs_edge_time_submat>max(t_vec))] = Inf
          }
          else{
            adjs_edge_time_submat = adjs_edge_time_submat + max(t_vec) - MAX_edge_time_submat
            adjs_edge_time_submat[which(adjs_edge_time_submat<0)] = Inf
          }
        }
      }
      
      cdf_array[q,l,] = ecdf(adjs_edge_time_submat)(t_vec)
      
      if(standardize & max(cdf_array[q,l,])>0){
        cdf_array[q,l,] = cdf_array[q,l,]/max(cdf_array[q,l,])
      }
    }
  }
  
  return(cdf_array)
}

