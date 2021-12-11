### Initialization of clusters and n0_mat randomly
get_init_v5 = function(edge_time_mat_list, N_clus, 
                       N_restart=1,
                       t_vec=seq(0, 200, length.out=1000))
{
  time_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  v_vec_list_best = list()
  n0_vec_list_best = list()
  membership_list_best = list()
  clusters_list_best = list()
  loss_best = Inf
  for (ind_restart in 1:N_restart) {
    ### Get init
    v_vec_list = list()
    n0_vec_list = list()
    membership_list = list()
    clusters_list = list()
    for (m in 1:N_subj) {
      
      # Initialize time shifts -------------------------------------------------------
      edge_time_mat = edge_time_mat_list[[m]]
      earliest_edge_time = apply(edge_time_mat, 1, function(row) min(row[which(row>1)]))
      
      n0_vec = (earliest_edge_time)/time_unit
      n0_vec[n0_vec==Inf] = 0
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
    
    ### Calculate loss
    if (N_restart==1) {
      loss = 0 ### When N_restart==1, no need to calculate loss
    } else {
      n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
      center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                 clusters_list = clusters_list, 
                                                 n0_mat_list = n0_mat_list, 
                                                 t_vec = t_vec)
      loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
                          n0_mat_list = n0_mat_list, 
                          clusters_list = clusters_list, 
                          center_cdf_array = center_cdf_array, 
                          t_vec = t_vec)$loss
    }
    
    ### Update init_best
    if(loss < loss_best){
      v_vec_list_best = v_vec_list
      n0_vec_list_best = n0_vec_list
      membership_list_best = membership_list
      clusters_list_best = clusters_list
      loss_best = loss
    }
  }
  

  return(list(membership_list=membership_list_best, 
              clusters_list=clusters_list_best, 
              n0_vec_list=n0_vec_list_best,
              v_vec_list=v_vec_list_best))
}