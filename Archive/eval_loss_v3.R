### Loss = Sum of |N_{m,i,k}-S^{v}F_{z_{m,i},k}|^2
### Consider shape(cdf) and conn_prob seperately
eval_loss_v3 = function(edge_time_mat_list, 
                        n0_mat_list, clusters_list, center_cdf_array,
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
    node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, 
                                            clusters = clusters, 
                                            n0_mat = n0_mat, 
                                            t_vec = t_vec)
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        for (i in clusters[[q]]) {
          center_cdf_1 = node_cdf_array[i,k,]
          center_cdf_2 = center_cdf_array[q,k,]
          conn_prob_1 = tail(center_cdf_1, 1)
          conn_prob_2 = tail(center_cdf_2, 1)
          center_cdf_1_normed = center_cdf_1 / conn_prob_1
          center_cdf_2_normed = center_cdf_2 / conn_prob_2
          
          dist_1 = sum((center_cdf_2_normed-center_cdf_1_normed)^2)
          dist_2 = (conn_prob_2-conn_prob_1)^2
          dist_tmp = dist_1 * (conn_prob_1^2+conn_prob_2^2)/2 + 
            dist_2 * ( sum((center_cdf_2_normed*I(center_cdf_2_normed<1))^2)+
                         sum((center_cdf_1_normed*I(center_cdf_1_normed<1))^2) )/2
          loss = loss + dist_tmp*length(clusters[[k]])
        }
        
      }
    }
    
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
