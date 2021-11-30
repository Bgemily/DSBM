
# obtain connecting pattern (cdf of shifted events) between each node and each cluster
get_node_cdf_array_v2 = function(edge_time_mat, clusters, 
                                 n0_mat=0, freq_trun = 15, 
                                 t_vec=seq(0, 50, 0.05)){  
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
      
      cdf_array[i,l, ] = tryCatch(ecdf(adjs_edge_time_submat)(t_vec),
                                  error=function(x)return(numeric(length(t_vec))))
      
      ### Smooth cdf_array by Fourier truncation
      cdf_ql = cdf_array[i,l,]
      ext_length = length(cdf_ql)%/%10
      cdf_ql_extend = c(rep(head(cdf_ql,1), ext_length),
                        cdf_ql,
                        rep(tail(cdf_ql,1), ext_length))
      fft_ql_extend = fft(cdf_ql_extend)
      fft_ql_extend_trun = c(head(fft_ql_extend, freq_trun+1),
                             rep(0, length(fft_ql_extend)-2*freq_trun-1),
                             tail(fft_ql_extend, freq_trun))
      cdf_ql_extend_trun = Re(fft(fft_ql_extend_trun, inverse = TRUE))
      cdf_ql_trun = cdf_ql_extend_trun[(1+ext_length):(length(cdf_ql)+ext_length)]
      cdf_array[i,l,] = cdf_ql_trun
    }
  }
  
  
  
  return(cdf_array)
}

