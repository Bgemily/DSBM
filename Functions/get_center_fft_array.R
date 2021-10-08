

### Obtain truncated fourier series (smoothed point process) for each pair of clusters 
### using multiple subjects
get_center_fft_array = function(edge_time_mat_list, clusters_list, freq_trun=5,
                                   n0_mat_list=NULL, t_vec=seq(0, 200, length.out=1000) )
{  
  time_unit = t_vec[2]-t_vec[1]
  N_subj = length(edge_time_mat_list)
  
  if (length(unique(sapply(clusters_list,length)))>1) {
    stop("The number of clusters is not the same across subjects.")
  }
  N_clus = length(clusters_list[[1]])
  
  center_fft_array = array(dim=c(N_clus, N_clus, 2*freq_trun+1))

  
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      ### Get shifted_edge_time_submat_list
      shifted_edge_time_submat_list = list()
      for (m in 1:N_subj) {
        edge_time_submat = edge_time_mat_list[[m]][clusters_list[[m]][[q]], clusters_list[[m]][[l]], drop=F]
        if (q==l) {
          diag(edge_time_submat) = NA
        }
        
        tau_submat = time_unit * n0_mat_list[[m]][clusters_list[[m]][[q]], clusters_list[[m]][[l]], drop=F]
        edge_time_submat = edge_time_submat - tau_submat
        
        shifted_edge_time_submat_list[[m]] = edge_time_submat
      }
      
      ### Get point_proc_fft_mat: Fourier series for N_ij's for i,j in clus q and l
      tmp = round(unlist(shifted_edge_time_submat_list)/time_unit)
      tmp = tmp[!is.na(tmp)] # Remove NA's due to t_{i,i}
      
      if (length(tmp)>0) {
        point_proc_fft_mat = matrix(0,nrow=length(tmp), ncol=2*freq_trun+1) 
        for (i in 1:nrow(point_proc_fft_mat)) {
          ### Turn edge time into a point process
          point_proc_vec = rep(0,length(t_vec))
          if (tmp[i]<Inf) {
            point_proc_vec[tmp[i]] = 1
          }
          ### Get normalized fourier series
          point_proc_vec_fft = fft(point_proc_vec) / length(t_vec)
          ### Truncate fourier series
          point_proc_fft_mat[i,] = c(tail(point_proc_vec_fft, freq_trun),
                                     head(point_proc_vec_fft, freq_trun+1))
        }
        
        
        center_fft_array[q,l,] = colMeans(point_proc_fft_mat)
      }
      else{
        center_fft_array[q,l,] = rep(0, 2*freq_trun+1) 
      }
      
      
    }
  }
  
  return(center_fft_array)
}

