
### Obtain truncated fourier series (smoothed shifted point process) between each node and each cluster 
### using multiple subjects

get_node_fft_array = function(edge_time_mat, clusters, 
                              n0_mat=edge_time_mat*0, freq_trun=5,
                              t_vec=seq(0, 200, length.out=1000))
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
      
      ###############################
      tmp = round(c(edge_time_submat)/time_unit)
      tmp = tmp[!is.na(tmp)]
      if (length(tmp)==0) {
        node_fft_array[i,k, ] = rep(0,2*freq_trun+1)
      }
      else{
        point_proc_fft_mat = matrix(0,nrow=length(tmp), ncol=2*freq_trun+1) 
        for (j in 1:nrow(point_proc_fft_mat)) {
          ### Turn edge time into a point process
          point_proc_vec = rep(0,length(t_vec))
          if (tmp[j]<Inf) {
            point_proc_vec[tmp[j]] = 1
          }
          ### Get normalized fourier series
          point_proc_vec_fft = fft(point_proc_vec) / length(t_vec)
          ### Truncate fourier series
          point_proc_fft_mat[j,] = c(tail(point_proc_vec_fft, freq_trun),
                                     head(point_proc_vec_fft, freq_trun+1))
        }
        
        node_fft_array[i,k, ] = colMeans(point_proc_fft_mat)
      }

    }
  }
  
  return(node_fft_array)
}

