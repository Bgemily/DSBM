### Loss = \sum_{i,j: t_{i,j}<\infty} T^{-1}*|N/N(T)-S^{v}F/F(T)|^2 + 0.1* \sum_{i,j: i\neq j} [N_{i,j}(T)-F_{z_i,z_j}(T)]^2

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
            if(edge_time_mat[i,j]<Inf & i!=j){
              loss_ij = mean((counting_proc_array[i,j,]/max(counting_proc_array[i,j,]+.Machine$double.eps) -
                                event_rate_array[i,j,]/max(event_rate_array[i,j,]+.Machine$double.eps))^2) +
                  0.1*(max(counting_proc_array[i,j,])-max(event_rate_array[i,j,]))^2
            } else if(edge_time_mat[i,j]==Inf & i!=j){
              loss_ij = 0.1*(max(counting_proc_array[i,j,])-max(event_rate_array[i,j,]))^2
            } else{
              loss_ij = 0
            }
            loss = loss + loss_ij
          }
        }
        
      }
    }
    
  }
  
  return(list(loss=loss))
}




