
### Obtain connecting pattern for each pair of clusters 
get_center_cdf_array_v2 = function(edge_time_mat_list, clusters_list, 
                                   n0_mat_list=NULL, freq_trun=Inf,
                                   t_vec=seq(0, 50, 0.05)){  
  time_unit = t_vec[2]-t_vec[1]
  N_subj = length(edge_time_mat_list)
  
  if (length(unique(sapply(clusters_list,length)))>1) {
    stop("The number of clusters is not the same across subjects.")
  }
  N_clus = length(clusters_list[[1]])
  
  ### Adjust edge times by time shifts
  adjs_edge_time_mat_list = list()
  for (m in 1:N_subj) {
    adjs_edge_time_mat = edge_time_mat_list[[m]] - time_unit*n0_mat_list[[m]]
    diag(adjs_edge_time_mat) = NA # remove (i,i)
    
    adjs_edge_time_mat_list[[m]] = adjs_edge_time_mat
  }
  
  
  cdf_array = array(dim=c(N_clus,N_clus,length(t_vec)))
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      if (sum(sapply(clusters_list, function(clusters)length(clusters[[q]])*length(clusters[[l]]))) <= 1) {
        cdf_array[q,l, ] = numeric(length(t_vec))
        next
      }
      
      adjs_edge_time_submat_list = list()
      for (m in 1:N_subj) {
        adjs_edge_time_submat = adjs_edge_time_mat_list[[m]][clusters_list[[m]][[q]], clusters_list[[m]][[l]], drop=F]
        adjs_edge_time_submat_list[[m]] = adjs_edge_time_submat
      }
      
      cdf_array[q,l,] = ecdf(unlist(adjs_edge_time_submat_list))(t_vec)
      
      ### Smooth cdf_array by Fourier truncation
      if (freq_trun<Inf){
        cdf_ql = cdf_array[q,l,]
        ext_length = length(cdf_ql)%/%10
        cdf_ql_extend = c(rep(head(cdf_ql,1), ext_length),
                          cdf_ql,
                          rep(tail(cdf_ql,1), ext_length))
        fft_ql_extend = fft(cdf_ql_extend) / length(cdf_ql_extend)
        fft_ql_extend_trun = c(head(fft_ql_extend, freq_trun+1),
                               rep(0, length(fft_ql_extend)-2*freq_trun-1),
                               tail(fft_ql_extend, freq_trun))
        cdf_ql_extend_trun = Re(fft(fft_ql_extend_trun, inverse = TRUE))
        cdf_ql_trun = cdf_ql_extend_trun[(1+ext_length):(length(cdf_ql)+ext_length)]
        cdf_array[q,l,] = cdf_ql_trun
      }
      
    }
  }
  
  return(cdf_array)
}

