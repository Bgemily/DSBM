### Initialization of clusters and n0_mat randomly
get_init_v5 = function(edge_time_mat_list, N_clus, 
                       N_restart=1,
                       t_vec=seq(0, 200, length.out=1000))
{
  time_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  ## TODO: Add restart and select best init by quadratic loss.
  
  v_vec_list = list()
  n0_vec_list = list()
  membership_list = list()
  clusters_list = list()
  for (m in 1:N_subj) {
    
    # Initialize time shifts -------------------------------------------------------
    edge_time_mat = edge_time_mat_list[[m]]
    earliest_edge_time = apply(edge_time_mat, 1, function(row) min(row[which(row>1)]))
    
    n0_vec = (earliest_edge_time)/time_unit
    n0_vec = runif(length(n0_vec), min=0, max=n0_vec)
    n0_vec = round(n0_vec)
    
    n0_vec_list[[m]] = n0_vec
    v_vec_list[[m]] = n0_vec*time_unit
    
 
    # Initialize clusters -----------------------------------------------------
    membership = sample(1:N_clus, size=N_node_vec[m], replace = TRUE)
    membership_list[[m]] = membership
    
    clusters = mem2clus(membership = membership, N_clus_min = N_clus)
    clusters_list[[m]] = clusters
    
  }
  
  
  
  
  return(list(membership_list=membership_list, 
              clusters_list=clusters_list, 
              n0_vec_list=n0_vec_list,
              v_vec_list=v_vec_list))
}