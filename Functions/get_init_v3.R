### Initialization of clusters and n0_mat
### Initialize time shifts by jittering true time shifts.
get_init_v3 = function(edge_time_mat_list, N_clus, v_true_list, 
                       jitter_time_rad = 0,
                       t_vec=seq(0, 200, length.out=1000))
{
  # print("####################")
  # print("get_init_v3: Initialize time shifts by jittering true time shifts.")
  # print("####################")
  
  time_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  

  n0_vec_list = list()
  membership_list = list()
  clusters_list = list()
  for (m in 1:N_subj) {
    
    # Initialize time shifts -------------------------------------------------------
    v_vec = v_true_list[[m]]
    
    ### Jitter time shifts
    v_vec = v_vec + stats::runif(length(v_vec), -jitter_time_rad, jitter_time_rad)
    if (sd(v_vec)!=0) {
      v_vec = v_vec / sd(v_vec) * sd(v_true_list[[m]]) ### Rescale
    }
    v_vec = v_vec - min(v_vec)
    
    n0_vec = v_vec / time_unit
    n0_vec = round(n0_vec)
    n0_vec = n0_vec - min(n0_vec)
    # n0_mat = n0_vec2mat(n0_vec = n0_vec)
    
    n0_vec_list[[m]] = n0_vec
    
    
    # Initialize clusters -----------------------------------------------------
    edge_time_mat = edge_time_mat_list[[m]]
    N_node = nrow(edge_time_mat)
    node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                           n0_mat = matrix(0,N_node,N_node), t_vec = t_vec)
    aligned_cdf_mat = t(sapply(1:N_node, function(i) shift_v2(f_origin=node_cdf_array[i,1,], 
                                                              n0=-n0_vec[i]) ))
    membership = cluster::pam(x=aligned_cdf_mat, k=N_clus, diss=FALSE, cluster.only=TRUE)
    membership_list[[m]] = membership

    clusters = mem2clus(membership)
    clusters_list[[m]] = clusters
    
  }
  
  
  
  
  return(list(membership_list=membership_list, 
              clusters_list=clusters_list, 
              n0_vec_list=n0_vec_list))
}

### Test
# edge_time_mat = matrix(1:4,2,2)
# edge_time_mat = kronecker(edge_time_mat, matrix(10,5,5))
# edge_time_mat_list = list(edge_time_mat, edge_time_mat)
# res = get_init_v2(edge_time_mat_list = edge_time_mat_list, N_clus = 2)
# res$clusters_list
# res$n0_vec_list
# res$membership_list
