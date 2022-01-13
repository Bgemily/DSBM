
### Update time shifts and connecting patterns alternatively

### Based on v7.1
### When getting rough estimates by aligning cdf's, iterate several times rather than only once.
### [WARNING: The gradient might be WRONG when there are multiple subjects!]

est_n0_vec_pdf = function(edge_time_mat_list, 
                         clusters_list, 
                         n0_vec_list=NULL, n0_mat_list=NULL,
                         freq_trun=5,
                         t_vec=seq(0,200,length.out=1000), step_size=200,
                         max_iter=50, epsilon=0.001, shrink=0.5, order_list=NULL, 
                         opt_radius=max(t_vec)/2,
                         ...)
{
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  N_clus = length(clusters_list[[1]])
  membership_list = lapply(clusters_list, clus2mem)
  
  if (is.null(n0_vec_list)) {
    res = get_init_v2(edge_time_mat_list=edge_time_mat_list, N_clus=N_clus, t_vec=t_vec)
    n0_vec_list = res$n0_vec_list
  }
  if (is.null(n0_mat_list)) {
    n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  }
  center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                          clusters_list = clusters_list, 
                                          freq_trun = freq_trun, 
                                          n0_mat_list = n0_mat_list, 
                                          t_vec = t_vec)
  
  
  # Update time shifts and connecting patterns alternatively ----------------
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  center_fft_array_update = center_fft_array_current = center_fft_array
  
  ### Set n0_min and n0_max
  n0_max_vec_list = n0_min_vec_list = list()
  for (m in 1:N_subj) {
    diag(edge_time_mat_list[[m]]) = Inf
    n0_max_vec = apply(edge_time_mat_list[[m]], 1, min)/t_unit
    n0_max_vec = pmin(n0_max_vec, n0_vec_list_current[[m]] + opt_radius/t_unit )
    n0_max_vec = round(n0_max_vec)
    
    
    n0_min_vec = 0*n0_max_vec
    n0_min_vec = pmax(n0_min_vec, n0_vec_list_current[[m]] - opt_radius/t_unit )
    n0_min_vec = round(n0_min_vec)
    
    
    n0_max_vec_list[[m]] = n0_max_vec
    n0_min_vec_list[[m]] = n0_min_vec
  }
  
  
  ### Refine estimation of time shifts using pdf
  n_iter = 1
  converge = FALSE
  while(!converge && n_iter<= max_iter)
  {
    
    ### Calculate the gradient at current time shifts
    gradient_vec_list = get_gradient(edge_time_mat_list = edge_time_mat_list, 
                                     clusters_list = clusters_list, 
                                     n0_vec_list = n0_vec_list_current, 
                                     n0_mat_list = n0_mat_list_current, 
                                     freq_trun = freq_trun, t_vec = t_vec)
    
    ### Update n0_vec_list_update
    for (m in 1:N_subj) {
      step_size_correct = step_size*(50/N_node_vec[1])
      n0_vec_list_update[[m]] = n0_vec_list_current[[m]] - step_size_correct*gradient_vec_list[[m]]  
      n0_vec_list_update[[m]] = round(n0_vec_list_update[[m]])
      
      ### Get n0_min and n0_max
      n0_max_vec = n0_max_vec_list[[m]]
      n0_min_vec = n0_min_vec_list[[m]]
      
      ### Force time shifts be between n0_min and n0_max
      for (i in 1:length(n0_vec_list_update[[m]])) {
        if(n0_vec_list_update[[m]][i] < n0_min_vec[i]){
          n0_vec_list_update[[m]][i] = n0_min_vec[i]
        }
        else if (n0_vec_list_update[[m]][i] > n0_max_vec[i]){
          n0_vec_list_update[[m]][i] = n0_max_vec[i]
        }
      }
      
      ### Set minimum shifted edge time to be zero
      n0_mat_tmp = n0_vec2mat(n0_vec = n0_vec_list_update[[m]])
      adjs_edge_time_mat = edge_time_mat_list[[m]] - n0_mat_tmp*t_unit
      global_timeshift = min(adjs_edge_time_mat)
      n0_vec_list_update[[m]] = n0_vec_list_update[[m]] + floor(abs(global_timeshift/t_unit))*sign(global_timeshift)

    }
    
    
    n0_mat_list_update = n0_vec2mat(n0_vec = n0_vec_list_update)
    ### Update connecting patterns 
    center_fft_array_update = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                                   clusters_list = clusters_list, 
                                                   freq_trun = freq_trun, 
                                                   n0_mat_list = n0_mat_list_update, 
                                                   t_vec = t_vec)
    
    ### Evaluate stopping criterion
    converge = sqrt(sum(abs(unlist(center_fft_array_update)-unlist(center_fft_array_current))^2))/
      sqrt(sum(abs(unlist(center_fft_array_current)+.Machine$double.eps)^2)) < epsilon
    
    
    ### *update -> *current
    n_iter = n_iter + 1
    n0_vec_list_current = n0_vec_list_update
    n0_mat_list_current = n0_mat_list_update
    center_fft_array_current = center_fft_array_update
  }
  
  if (n_iter>max_iter) {
    message("[est_n0_vec_pdf]: Reached maximum iteration number ", max_iter)
  }
  
  
  ### Compute time shift matrix and connecting patterns
  n0_vec_list = n0_vec_list_current
  n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  
  v_vec_list = lapply(n0_vec_list, function(n0_vec)n0_vec*t_unit)
  v_mat_list = lapply(n0_mat_list, function(n0_mat)n0_mat*t_unit)
  
  center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                          clusters_list = clusters_list, 
                                          freq_trun = freq_trun, 
                                          n0_mat_list = n0_mat_list, 
                                          t_vec = t_vec)
  
  return(list(center_fft_array=center_fft_array, 
              n0_vec_list=n0_vec_list, n0_mat_list=n0_mat_list,
              v_vec_list=v_vec_list, v_mat_list=v_mat_list))
}




# Test --------------------------------------------------------------------

# res = generate_network2_v3(N_subj = 1,N_node_vec = c(180), total_time = 200,conn_patt_var = 1,
#                            conn_patt_sep = 1.3, conn_prob_mean = 1, conn_prob_rad = 0,
#                            time_shift_mean_vec = rep(5,3),time_shift_rad = 5,)
# edge_time_mat_list = res$edge_time_mat_list
# clusters_list = res$clus_true_list
# time_shift_list = res$time_shift_list
# pdf_true_array = res$pdf_true_array
# plot_pdf_array_v2(pdf_true_array)
# 
# 
# truth = (time_shift_list[[1]])
# 
# edge_time_mat_list = list(edge_time_mat_list[[1]][1:60,1:60])
# clusters_list[[1]][2:3] = NULL
# truth = truth[1:60]
# 
# tmp = est_n0_vec_v7(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list,
#                     freq_trun = 10,step_size = 5, max_iter = 20,
#                     t_vec = seq(0,200,1),
#                     )
#                     # n0_vec_list = list(truth/0.2+runif(length(truth),0,20)))
# plot(truth, tmp$v_vec_list[[1]])
# abline(a=0,b=1,col=2)
# 
# 
# 
# tmp2 = est_n0_vec_cdf(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list,
#                        t_vec = seq(0,200,1),
#                        )
#                       # n0_vec_list = list(truth/0.2+runif(length(truth),0,20)), 
#                       # n0_mat_list = n0_vec2mat(list(truth/0.2+runif(length(truth),0,20))) )
# plot(truth, tmp2$v_vec_list[[1]])
# abline(a=0,b=1,col=2)

