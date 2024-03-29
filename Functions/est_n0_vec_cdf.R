
### Centering step: Update time shifts and connecting patterns alternatively
### Force n0 to be less than earliest edge time.

est_n0_vec_cdf = function(edge_time_mat_list, 
                         clusters_list, center_cdf_array=NULL, 
                         n0_vec_list=NULL, n0_mat_list=NULL,
                         freq_trun = Inf,
                         t_vec=seq(0,200,length.out=1000), step_size=0.02,
                         max_iter=5, epsilon=0.001)
{
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  N_clus = length(clusters_list[[1]])

  if (is.null(n0_vec_list)) {
    res = get_init_v2(edge_time_mat_list=edge_time_mat_list, N_clus=N_clus, t_vec=t_vec)
    n0_vec_list = res$n0_vec_list
    n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  }
  
  if (is.null(center_cdf_array)) {
    center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                               clusters_list = clusters_list, 
                                               n0_mat_list = n0_mat_list, 
                                               freq_trun = freq_trun, 
                                               t_vec = t_vec)
  }
  

  # Update time shifts and connecting patterns alternatively ----------------
  n_iter = 1
  converge = FALSE
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  center_cdf_array_update = center_cdf_array_current = center_cdf_array
  while(!converge && n_iter<= max_iter)
  {
    ### Update time shifts given connecting patterns and order of time shifts
    for (m in 1:N_subj) {
      n0_vec_tmp = n0_vec_list_current[[m]]
      N_node = N_node_vec[[m]]
      
      ### Get order matrix with (i,j)-th entry being 1 if v_j<=v_i and NA otherwise
      n0_vec_tmp[which.min(n0_vec_tmp)] = min(n0_vec_tmp[-which.min(n0_vec_tmp)])
      order_mat_current = matrix(n0_vec_tmp, nrow=N_node, ncol=N_node) - 
        matrix(n0_vec_tmp, nrow=N_node, ncol=N_node, byrow=TRUE) >= 0
      order_mat_current = order_mat_current*1 ### Convert logical entries to numeric
      diag(order_mat_current) = 0
      order_mat_current[order_mat_current==0] = NA
      
      ### Get aggregated counting processes \bar N_{m,i,k}
      edge_time_mat_tmp = edge_time_mat_list[[m]] * order_mat_current
      clusters_tmp = clusters_list[[m]]
      n0_mat_tmp = n0_mat_list_current[[m]]
      node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat_tmp, 
                                             clusters = clusters_tmp, 
                                             n0_mat = n0_mat_tmp*0, 
                                             freq_trun = freq_trun, 
                                             t_vec = t_vec)
      
      ### Get weight matrix N_node*N_clus, (i,l)-th entry is proportional to |{j: v_j<=v_i && z_j==l}|
      order_mat_current[is.na(order_mat_current)] = 0
      exist_edge_mat = edge_time_mat_list[[m]]<Inf
      diag(exist_edge_mat) = 0
      weight_mat = matrix(nrow=N_node, ncol=N_clus)
      for (q in 1:N_clus) {
        weight_mat[,q] = rowSums(as.matrix((order_mat_current*exist_edge_mat)[,clusters_tmp[[q]]]))
      }
      
      ### Update time shifts
      for (q in 1:N_clus) {
        for (i in clusters_tmp[[q]]) {
          f_target_list = lapply(1:N_clus, function(l) node_cdf_array[i,l, ])
          f_origin_list = lapply(1:N_clus, function(l) center_cdf_array_current[q,l, ])
          
          ##### Eliminate the difference between the constant values at the right end of curves.
          f_origin_list = lapply(1:length(clusters_tmp), function(l) f_origin_list[[l]]/max(max(f_origin_list[[l]]), 1e-6))
          f_target_list = lapply(1:length(clusters_tmp), function(l) f_target_list[[l]]/max(max(f_target_list[[l]]), 1e-6))
          
          weights = weight_mat[i,]
          if(sum(weights)==0)
            weights = NULL
          else
            weights = weights/sum(weights)
          
          
          n0_max = min(edge_time_mat_list[[m]][i,])/t_unit
          n0_max = round(n0_max)
          n0_max = min(n0_max, length(t_vec)-1)
          n0_vec_tmp[i] = align_multi_curves_gd_v2(f_origin_list = f_origin_list, f_target_list = f_target_list,
                                                   step_size = step_size,
                                                   t_unit = t_unit, weights = weights,
                                                   n0_max = n0_max)$n0
          
        }
      }
      n0_vec_list_update[[m]] = n0_vec_tmp
      
      ### Set minimum shifted edge time to be one
      n0_mat_tmp = n0_vec2mat(n0_vec = n0_vec_list_update[[m]])
      adjs_edge_time_mat = edge_time_mat_list[[m]] - n0_mat_tmp*t_unit
      global_timeshift = min(adjs_edge_time_mat)-1
      n0_vec_list_update[[m]] = n0_vec_list_update[[m]] + floor(abs(global_timeshift/t_unit))*sign(global_timeshift)
      
    }

    ### Update connecting patterns 
    n0_mat_list_update = n0_vec2mat(n0_vec = n0_vec_list_update)
    
    center_cdf_array_update = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                      clusters_list = clusters_list, 
                                                      n0_mat_list = n0_mat_list_update, 
                                                      freq_trun = freq_trun, 
                                                      t_vec = t_vec)
    ### Evaluate stopping criterion
    converge = sqrt(sum((unlist(center_cdf_array_update)-unlist(center_cdf_array_current))^2))/
      sqrt(sum((unlist(center_cdf_array_current)+.Machine$double.eps)^2)) < epsilon
    
    ### *update -> *current
    n_iter = n_iter + 1
    n0_vec_list_current = n0_vec_list_update
    n0_mat_list_current = n0_mat_list_update
    center_cdf_array_current = center_cdf_array_update
  }
  
  ### Compute time shift matrix and connecting patterns
  n0_vec_list = n0_vec_list_current
  n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  
  v_vec_list = lapply(n0_vec_list, function(n0_vec)n0_vec*t_unit)
  v_mat_list = lapply(n0_mat_list, function(n0_mat)n0_mat*t_unit)
  
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, 
                                             freq_trun = Inf, 
                                             t_vec = t_vec)
  
  return(list(center_cdf_array=center_cdf_array, 
              n0_vec_list=n0_vec_list, n0_mat_list=n0_mat_list,
              v_vec_list=v_vec_list, v_mat_list=v_mat_list))
}


