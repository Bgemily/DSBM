
### Obtain truncated fourier series (smoothed shifted point process) between each node and each cluster 
### using multiple subjects

get_node_fft_array = function(edge_time_mat, clusters, 
                              n0_mat=edge_time_mat*0, freq_trun=15,
                              t_vec=seq(0, 200, length.out=1000),
                              rmv_conn_prob=FALSE)
{  
  time_unit = t_vec[2]-t_vec[1]
  N_node = nrow(edge_time_mat)
  N_clus = length(clusters)
  
  ### Set the edge time t_{i,i}'s to be NA
  diag(edge_time_mat) = NA
  
  node_fft_array = array(dim=c(N_node,N_clus,2*freq_trun+1))
  for (i in 1:N_node) {
    for (k in 1:N_clus) {
      edge_time_submat = edge_time_mat[i, clusters[[k]], drop=F]
      
      tau_submat = time_unit * n0_mat[i, clusters[[k]], drop=F]
      edge_time_submat = edge_time_submat - tau_submat
      
      #### Get node_fft_array[i,k, ]
      tmp = c(edge_time_submat)
      tmp = tmp[!is.na(tmp)] # Remove NA's due to t_{i,i}
      if(rmv_conn_prob){
        tmp = tmp[tmp<Inf] ### Remove non-exist edges
      }
      
      if (length(tmp)==0) {
        node_fft_array[i,k, ] = rep(0,2*freq_trun+1)
      }
      else{
        if (min(tmp)<0) {
          tmp = tmp - min(tmp)
        }
        fft_res = get_adaptive_fft(event_time_vec = tmp, 
                                   freq_trun_max = freq_trun, 
                                   t_vec = t_vec)
        node_fft_array[i,k, ] = fft_res$fft_vec_best
      }

    }
  }
  
  return(node_fft_array)
}

