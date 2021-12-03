
### Update time shifts and connecting patterns alternatively

### Based on v4
### Force n0 to be less than earliest edge time.
### Normalize node_cdf_array when estimating time shifts

est_n0_vec_v4.1 = function(edge_time_mat_list, 
                         clusters_list, center_cdf_array=NULL, 
                         n0_vec_list=NULL, n0_mat_list=NULL,
                         freq_trun = Inf,
                         t_vec=seq(0,200,length.out=1000), step_size=0.02,
                         max_iter=5, epsilon=0.001, order_list=NULL)
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
  
  # browser(); mean(center_cdf_array)
  
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
      weight_mat = matrix(nrow=N_node, ncol=N_clus)
      for (q in 1:N_clus) {
        weight_mat[,q] = rowSums(as.matrix(order_mat_current[,clusters_tmp[[q]]]))
      }
      
      ### Update time shifts
      for (q in 1:N_clus) {
        for (i in clusters_tmp[[q]]) {
          f_target_list = lapply(1:N_clus, function(l) node_cdf_array[i,l, ])
          f_origin_list = lapply(1:N_clus, function(l) center_cdf_array_current[q,l, ])
          
          ##### V1: Eliminate the difference between the constant values at the right end of curves.
          f_origin_list = lapply(1:length(clusters_tmp), function(l) f_origin_list[[l]]/max(max(f_origin_list[[l]]), 1e-6))
          f_target_list = lapply(1:length(clusters_tmp), function(l) f_target_list[[l]]/max(max(f_target_list[[l]]), 1e-6))
          
          ##### V2
          # blank
          ####################
          
          weights = weight_mat[i,]
          if(sum(weights)==0)
            weights = NULL
          else
            weights = weights/sum(weights)
          
          
          if (!is.null(order_list)) {
            position = which(order_list[[m]]==i)
            n0_min = ifelse(test = position>=2, 
                            yes = n0_vec_tmp[order_list[[m]][position-1]],
                            no = 0)
            n0_max = ifelse(test = position<=(N_node-1), 
                            yes = n0_vec_tmp[order_list[[m]][position+1]],
                            no = length(t_vec))
            n0_vec_tmp[i] = align_multi_curves_gd_v2(f_origin_list = f_origin_list, f_target_list = f_target_list,
                                                     step_size = step_size,
                                                     t_unit = t_unit, weights = weights,
                                                     n0_min = n0_min, n0_max = n0_max)$n0
          }
          else{
            n0_max = min(edge_time_mat_list[[m]][i,])/t_unit
            n0_max = round(n0_max)
            n0_max = min(n0_max, length(t_vec)-1)
            n0_vec_tmp[i] = align_multi_curves_gd_v2(f_origin_list = f_origin_list, f_target_list = f_target_list,
                                                     step_size = step_size,
                                                     t_unit = t_unit, weights = weights,
                                                     n0_max = n0_max)$n0
            
          }
          
        }
      }
      n0_vec_list_update[[m]] = n0_vec_tmp
      
      ### Set minimum shifted edge time to be zero
      n0_mat_tmp = n0_vec2mat(n0_vec = n0_vec_list_update[[m]])
      adjs_edge_time_mat = edge_time_mat_list[[m]] - n0_mat_tmp*t_unit
      global_timeshift = min(adjs_edge_time_mat)
      n0_vec_list_update[[m]] = n0_vec_list_update[[m]] + floor(abs(global_timeshift/t_unit))*sign(global_timeshift)
      
    }
    ### Debug
    # browser()
    # plot(n0_vec_list_current[[1]], n0_vec_list_update[[1]])
    
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
  
  if (n_iter>max_iter) {
    message("[est_n0_vec_v4.1]: Reached maximum iteration number ", max_iter)
  }
  
  ### Compute time shift matrix and connecting patterns
  n0_vec_list = n0_vec_list_current
  n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  
  v_vec_list = lapply(n0_vec_list, function(n0_vec)n0_vec*t_unit)
  v_mat_list = lapply(n0_mat_list, function(n0_mat)n0_mat*t_unit)
  
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, 
                                             freq_trun = freq_trun, 
                                             t_vec = t_vec)
  
  return(list(center_cdf_array=center_cdf_array, 
              n0_vec_list=n0_vec_list, n0_mat_list=n0_mat_list,
              v_vec_list=v_vec_list, v_mat_list=v_mat_list))
}




# Test --------------------------------------------------------------------

# res = generate_network2_v2(N_subj = 1,N_node_vec = c(30), N_clus = 3,total_time = 200,conn_patt_var = 1,
#                            conn_patt_sep = 1.5, conn_prob_mean = 0.8,conn_prob_rad = 0,
#                            time_shift_mean_vec = rep(10,3),time_shift_rad = 10,)
# 
# edge_time_mat_list = res$edge_time_mat_list
# cdf_true_array = res$cdf_true_array
# pdf_true_array = res$pdf_true_array
# clusters_list = res$clus_true_list
# time_shift_list = res$time_shift_list
# truth = (time_shift_list[[1]])
# 
# tmp = est_n0_vec_v4(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list)
# tmp$v_vec_list
# plot(truth, tmp$v_vec_list[[1]])
# 
# est = est_n0_vec(edge_time_mat = edge_time_mat_list[[1]], clusters = clusters_list[[1]], t_vec = seq(0,200,length.out=1000))
# plot(truth, est)
