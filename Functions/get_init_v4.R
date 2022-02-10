### Initialize of cluster memberships and time shifts
### Initialize time shifts by earliest edge time.
get_init_v4 = function(edge_time_mat_list, N_clus, 
                       freq_trun=Inf,
                       t_vec=seq(0, 200, length.out=1000))
{

  time_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  v_vec_list = list()
  n0_vec_list = list()
  membership_list = list()
  clusters_list = list()
  for (m in 1:N_subj) {
    
    # Initialize time shifts -------------------------------------------------------
    
    edge_time_mat = edge_time_mat_list[[m]]
    earliest_edge_time = apply(edge_time_mat, 1, function(row) min(row[which(row>1)]))
    
    n0_vec = (earliest_edge_time)/time_unit
    n0_vec = round(n0_vec)
    n0_vec[n0_vec==Inf] = 0
    
    n0_vec = n0_vec - min(n0_vec)
    n0_vec_list[[m]] = n0_vec
    v_vec_list[[m]] = n0_vec*time_unit
    
    
    
    # Initialize clusters -----------------------------------------------------
    edge_time_mat = edge_time_mat_list[[m]]
    N_node = nrow(edge_time_mat)
    node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, 
                                           clusters = list(c(1:N_node)), 
                                           n0_mat = matrix(n0_vec,nrow=N_node,ncol=N_node), 
                                           freq_trun = freq_trun,
                                           t_vec = t_vec)
    aligned_cdf_mat = t(sapply(1:N_node, function(i) node_cdf_array[i,1,]/max(node_cdf_array[i,1,]+.Machine$double.eps) ))
    membership = cluster::pam(x=aligned_cdf_mat, k=N_clus, diss=FALSE, cluster.only=TRUE)
    membership_list[[m]] = membership

    clusters = mem2clus(membership = membership, N_clus_min = N_clus)
    clusters_list[[m]] = clusters
    
  }
  
  
  
  
  return(list(membership_list=membership_list, 
              clusters_list=clusters_list, 
              n0_vec_list=n0_vec_list,
              v_vec_list=v_vec_list))
}

