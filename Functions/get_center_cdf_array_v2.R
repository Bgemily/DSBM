
### Obtain connecting pattern for each pair of clusters using multiple subjects
get_center_cdf_array_v2 = function(edge_time_mat_list, clusters_list, 
                                   n0_mat_list=NULL, freq_trun=15,
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
      cdf_ql = cdf_array[q,l,]
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
      cdf_array[q,l,] = cdf_ql_trun
    }
  }
  
  return(cdf_array)
}


# ### Test
# edge_time_mat1 = kronecker(matrix(1:4,2,2), matrix(10,5,5))
# edge_time_mat2 = kronecker(matrix(1:4,2,2), matrix(10,3,3))+1
# clusters1 = list(1:5,6:10)
# clusters2 = list(1:3,4:6)
# edge_time_mat_list = list(edge_time_mat1)
# clusters_list = list(clusters1)
# n0_mat_list = list(edge_time_mat1*0+1/0.05)
# 
# res = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list,
#                               clusters_list = clusters_list,
#                               freq_trun = 20,
#                               t_vec = seq(0,50,length.out=200),
#                               n0_mat_list = n0_mat_list)
# par(mfrow=c(2,2))
# for (q in 1:2) {
#   for (k in 1:2) {
#     plot(res[q,k,], type='l')
#   }
# }
# par(mfrow=c(1,1))
