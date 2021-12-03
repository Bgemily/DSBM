### Loss = Sum of |N-S^{v}F|^2

eval_loss_v2 = function(edge_time_mat_list, 
                        n0_mat_list, clusters_list, center_cdf_array,
                        freq_trun = Inf,
                        t_vec=seq(0,200,length.out=1000)){

  N_subj = length(edge_time_mat_list)
  N_clus = length(clusters_list[[1]])
  t_unit = t_vec[2]-t_vec[1]
  
  
  ### Evaluate loss (up to a constant)
  loss = 0
  for (m in 1:N_subj) {
    edge_time_mat = edge_time_mat_list[[m]]
    clusters = clusters_list[[m]]
    n0_mat = n0_mat_list[[m]]
    
    N_node = nrow(edge_time_mat)
    counting_proc_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, 
                                                clusters = as.list(1:N_node), 
                                                n0_mat = matrix(0,nrow=N_node,ncol=N_node),
                                                freq_trun = freq_trun,
                                                t_vec = t_vec)
    event_rate_array = array(dim = c(N_node,N_node,length(t_vec)))
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        
        for (i in clusters[[q]]) {
          for (j in clusters[[k]]) {
            event_rate_array[i,j,] = shift_v2(f_origin = center_cdf_array[q,k,], n0 = n0_mat[i,j])
            if (i==j) {
              event_rate_array[i,j,] = counting_proc_array[i,j,] ### To make sure (N-S^{v}F)_{ii}==0
            }
          }
        }
        
      }
    }
    
    loss = loss + sum((counting_proc_array-event_rate_array)^2)
    
  }
  
  return(list(loss=loss))
}





# Test --------------------------------------------------------------------
# 
# edge_time_mat = kronecker(matrix(1,2,2), matrix(100,5,5))
# t_vec = seq(0,200,0.2)
# t_unit = t_vec[2] - t_vec[1]
# n0_mat = (edge_time_mat*0+5)/t_unit
# edge_time_mat = edge_time_mat + 5
# cdf_true_array = ecdf(100)(t_vec)
# clusters = list(1:5,6:10)
# 
# edge_time_mat_list = list(edge_time_mat,edge_time_mat)
# n0_mat_list = list(n0_mat, n0_mat)
# clusters_list = list(clusters,clusters)
# center_cdf_array = array(dim=c(2,2,length(t_vec)))
# center_cdf_array[1,1,] = center_cdf_array[2,2,] = center_cdf_array[1,2,] = center_cdf_array[2,1,] = cdf_true_array
# 
# 
# loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list,
#                    n0_mat_list = n0_mat_list,
#                    clusters_list = clusters_list,
#                    center_cdf_array = center_cdf_array, t_vec = t_vec)$loss
# 
# 
# 
