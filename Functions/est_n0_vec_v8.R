
### Update time shifts and connecting patterns alternatively

### Based on v4.1
### Switch from cdf to pdf: pairwise alignment
###[TO DO: make sure there is no cdf or pdf in this file]

est_n0_vec_v8 = function(edge_time_mat_list, 
                         clusters_list, center_fft_array=NULL, 
                         n0_vec_list=NULL, n0_mat_list=NULL,
                         freq_trun=5,
                         t_vec=seq(0,200,length.out=1000), step_size=0.02,
                         max_iter=5, epsilon=0.001, order_list=NULL, 
                         opt_radius=max(t_vec)/2,...)
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
  
  if (is.null(center_fft_array)) {
    center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                            clusters_list = clusters_list, 
                                            freq_trun = freq_trun, 
                                            n0_mat_list = n0_mat_list, 
                                            t_vec = t_vec)
  }
  
  
  # Update time shifts and connecting patterns alternatively ----------------
  n_iter = 1
  converge = FALSE
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  center_fft_array_update = center_fft_array_current = center_fft_array
  while(!converge && n_iter<= max_iter)
  {
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
      node_fft_array = get_node_fft_array(edge_time_mat = edge_time_mat_tmp, 
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
          fft_target_list = lapply(1:N_clus, function(l) node_fft_array[i,l, ])
          fft_origin_list = lapply(1:N_clus, function(l) center_fft_array_current[q,l, ])
          
          ### fft_target_list: tail-fft(freq_trun) + head-fft(freq_trun+1)
          ### Need: head-fft(freq_trun+1) + rep(0,x) + tail-fft(freq_trun)
          fft_target_list = lapply(fft_target_list, 
                                   function(fft_vec)c(tail(fft_vec,freq_trun+1), 
                                                      rep(0,length(t_vec)-2*freq_trun-1),
                                                      head(fft_vec,freq_trun)))
          fft_origin_list = lapply(fft_origin_list, 
                                   function(fft_vec)c(tail(fft_vec,freq_trun+1), 
                                                      rep(0,length(t_vec)-2*freq_trun-1),
                                                      head(fft_vec,freq_trun)))
          
          
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
            stop("[est_n0_vec_v8:] The function align_multi_fft_gd might have bugs.")
            n0_vec_tmp[i] = align_multi_fft_gd(fft_origin_list = fft_origin_list, fft_target_list = fft_target_list,
                                                     step_size = step_size,
                                                     freq_trun = freq_trun,
                                                     weights = weights, t_vec=t_vec,
                                                     n0_min = n0_min, n0_max = n0_max)$n0
          }
          else{
            n0_max = n0_max_vec_list[[m]][i]
            n0_min = n0_min_vec_list[[m]][i]
            
            f_origin_list = lapply(fft_origin_list,function(fft)Re(fft(fft, inverse=TRUE)) )
            f_target_list = lapply(fft_target_list,function(fft)Re(fft(fft, inverse=TRUE)) )
            n0_vec_tmp[i] = align_multi_curves_gd_v2(f_origin_list = f_origin_list,
                                                     f_target_list = f_target_list,
                                                     step_size = step_size,
                                                     weights = weights,
                                                     t_unit = t_vec[2]-t_vec[1],
                                                     n0 = n0_vec_list_current[[m]][i],
                                                     n0_min = 0, n0_max = n0_max)$n0
          }
          
        }
      }
      
      n0_vec_tmp = n0_vec_tmp - min(n0_vec_tmp)
      n0_vec_list_update[[m]] = n0_vec_tmp
      
    }
    ### Debug
    # browser()
    # plot(n0_vec_list_current[[1]], n0_vec_list_update[[1]])
    
    ### Update connecting patterns 
    n0_mat_list_update = n0_vec2mat(n0_vec = n0_vec_list_update)
    
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
# tmp2 = est_n0_vec_v4.1(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list,
#                        t_vec = seq(0,200,1),
#                        )
#                       # n0_vec_list = list(truth/0.2+runif(length(truth),0,20)), 
#                       # n0_mat_list = n0_vec2mat(list(truth/0.2+runif(length(truth),0,20))) )
# plot(truth, tmp2$v_vec_list[[1]])
# abline(a=0,b=1,col=2)

