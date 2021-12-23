
### Loss = Sum of |dN/dt-S^{v}f|^2 in the frequency domain using cut-off frequency
### Based on v2
### Switch from cdf to pdf, switch from time domain to frequency domain

eval_loss_pdf = function(edge_time_mat_list, clusters_list, 
                        n0_vec_list=NULL, n0_mat_list=NULL,
                        center_fft_array=NULL,
                        freq_trun = 5, t_vec=seq(0,200,length.out=1000)){

  N_subj = length(edge_time_mat_list)
  N_clus = length(clusters_list[[1]])
  t_unit = t_vec[2]-t_vec[1]
  
  if (is.null(n0_mat_list)) {
    n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  }
  
  if (is.null(center_fft_array)) {
    center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                            clusters_list = clusters_list, 
                                            n0_mat_list = n0_mat_list, 
                                            freq_trun = freq_trun, 
                                            t_vec = t_vec)
  }
  
  
  ### Evaluate loss (up to a constant)
  loss = 0
  for (m in 1:N_subj) {
    edge_time_mat = edge_time_mat_list[[m]]
    clusters = clusters_list[[m]]
    n0_mat = n0_mat_list[[m]]
    
    N_node = nrow(edge_time_mat)
    
    ### Calculate the Fourier series for shifted point processes, i.e. S^\max(v_i,v_j) N_{i,j}'s
    point_proc_fft_array = get_node_fft_array(edge_time_mat = edge_time_mat, 
                                            clusters = as.list(1:N_node), 
                                            n0_mat = n0_mat, t_vec = t_vec, 
                                            freq_trun = freq_trun)
    intensity_fft_array = array(dim = c(N_node,N_node,2*freq_trun+1))
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        
        for (i in clusters[[q]]) {
          for (j in clusters[[k]]) {
            intensity_fft_array[i,j,] = center_fft_array[q,k,]
            if (i==j) {
              intensity_fft_array[i,j,] = point_proc_fft_array[i,j,] = 0 ### To make sure (N-S^{v}F)_{ii}==0
            }
          }
        }
        
      }
    }
    
    loss = loss + sum((abs(point_proc_fft_array-intensity_fft_array))^2) * length(t_vec)
    
  }
  
  return(list(loss=loss))
}





# Test --------------------------------------------------------------------

# edge_time_mat = kronecker(matrix(1,2,2), matrix(100,5,5))
# t_vec = seq(0,200,0.2)
# t_unit = t_vec[2] - t_vec[1]
# n0_mat = (edge_time_mat*0+5)/t_unit
# edge_time_mat = edge_time_mat + 5
# 
# freq_trun = 5
# tmp = 0*t_vec
# tmp[500] = 1
# fft_true_array = fft(tmp)/length(t_vec)
# fft_true_array = c(tail(fft_true_array, freq_trun),
#                    head(fft_true_array, freq_trun+1))
# 
# clusters = list(1:5,6:10)
# 
# edge_time_mat_list = list(edge_time_mat,edge_time_mat)
# n0_mat_list = list(n0_mat, n0_mat)
# clusters_list = list(clusters,clusters)
# center_fft_array = array(dim=c(2,2,2*freq_trun+1))
# center_fft_array[1,1,] = center_fft_array[2,2,] = center_fft_array[1,2,] = center_fft_array[2,1,] = fft_true_array
# 
# 
# eval_loss_pdf(edge_time_mat_list = edge_time_mat_list,
#              n0_mat_list = lapply(n0_mat_list, function(x)x+rnorm(length(x),sd = 10)),
#              clusters_list = clusters_list,
#              center_fft_array = center_fft_array,
#              t_vec = t_vec)$loss



