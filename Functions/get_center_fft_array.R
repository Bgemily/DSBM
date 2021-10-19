

### Obtain truncated fourier series (smoothed point process) for each pair of clusters 
### using multiple subjects
get_center_fft_array = function(edge_time_mat_list, clusters_list, 
                                n0_mat_list=NULL, freq_trun=15,
                                t_vec=seq(0, 200, length.out=1000) )
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
      
      ### Get center_fft_array[q,l,]
      tmp = unlist(shifted_edge_time_submat_list)
      tmp = tmp[!is.na(tmp)] # Remove NA's due to t_{i,i}
      
      if (length(tmp)>0) {
        if(min(tmp)<0)
          tmp = tmp-min(tmp)
        fft_res = get_adaptive_fft(event_time_vec = tmp, 
                                   freq_trun_max = freq_trun, 
                                   t_vec = t_vec)
        center_fft_array[q,l,] = fft_res$fft_vec_best
      }
      else{
        center_fft_array[q,l,] = rep(0, 2*freq_trun+1) 
      }
      
      
    }
  }
  
  return(center_fft_array)
}

