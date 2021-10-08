
### Update time shifts and connecting patterns alternatively

### Based on v7
### When aligning pdfs, optimization region's radius is controlled by a tuning parameter. 
### [WARNING: The gradient might be WRONG when there are multiple subjects!]

est_n0_vec_v7.1 = function(edge_time_mat_list, 
                         clusters_list, 
                         n0_vec_list=NULL, n0_mat_list=NULL,
                         freq_trun=5,
                         t_vec=seq(0,200,length.out=1000), step_size=NULL,
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

  
  # Update time shifts and connecting patterns alternatively ----------------
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  center_cdf_array_current = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                       clusters_list = clusters_list, 
                                                       n0_mat_list = n0_mat_list, t_vec = t_vec)
  center_cdf_array_update = center_cdf_array_current
  

  ### Get rough estimation of time shifts using cdf
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
    node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat_tmp, clusters = clusters_tmp,
                                           n0_mat = n0_mat_tmp*0, t_vec = t_vec)

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
                                                   step_size = 0.02,
                                                   t_unit = t_unit, weights = weights,
                                                   n0_min = n0_min, n0_max = n0_max)$n0
        }
        else{
          n0_max = min(edge_time_mat_list[[m]][i,])/t_unit
          n0_max = round(n0_max)
          n0_max = min(n0_max, length(t_vec)-1)
          n0_vec_tmp[i] = align_multi_curves_gd_v2(f_origin_list = f_origin_list, f_target_list = f_target_list,
                                                   step_size = 0.02,
                                                   t_unit = t_unit, weights = weights,
                                                   n0_max = n0_max)$n0

        }

      }
    }

    # n0_vec_tmp = n0_vec_tmp - min(n0_vec_tmp)
    n0_vec_list_update[[m]] = n0_vec_tmp

  }

  ### *update -> *current
  n0_vec_list_current = n0_vec_list_update
  n0_mat_list_current = n0_mat_list_update

  
  ### Evaluate loss 
  loss_history = c()
  loss = eval_loss_v4(edge_time_mat_list = edge_time_mat_list,
                      n0_vec_list = n0_vec_list_current,
                      clusters_list = clusters_list,
                      freq_trun = freq_trun,
                      t_vec = t_vec)$loss
  loss_history = c(loss_history, loss)
  
  
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
  
  
  ### Set initial learning rate
  gradient_vec_list = get_gradient(edge_time_mat_list = edge_time_mat_list, 
                                   clusters_list = clusters_list, 
                                   n0_vec_list = n0_vec_list_current, 
                                   n0_mat_list = n0_mat_list_current, 
                                   freq_trun = freq_trun, t_vec = t_vec)
  step_size = length(t_vec) / mean(sapply(gradient_vec_list, function(vec) sqrt(sum(vec^2)) ))
  
  
  ### Refine estimation of time shifts using pdf
  n_iter = 1
  converge = FALSE
  loss_current = loss
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
      n0_vec_list_update[[m]] = n0_vec_list_current[[m]] - step_size*gradient_vec_list[[m]]  
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
      

      ### Set the minimum time shift as zero
      n0_vec_list_update[[m]] = n0_vec_list_update[[m]] - min(n0_vec_list_update[[m]])
      
    }
    
    
    ### Evaluate loss 
    loss = eval_loss_v4(edge_time_mat_list = edge_time_mat_list,
                        n0_vec_list = n0_vec_list_update,
                        clusters_list = clusters_list,
                        freq_trun = freq_trun,
                        t_vec = t_vec)$loss
    loss_update = loss
    
    
    ### Adjust step size
    N_shrink = 0
    while (loss_update>loss_current & N_shrink<=10) {
      ### Shrink step_size
      N_shrink = N_shrink+1
      step_size = step_size*shrink
      
      ### Estimate n0_vec_list_update based on shrunk step size
      for (m in 1:N_subj) {
        n0_vec_list_update[[m]] = n0_vec_list_current[[m]] - step_size*gradient_vec_list[[m]]  
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
        
        ### Set the minimum time shift as zero
        n0_vec_list_update[[m]] = n0_vec_list_update[[m]] - min(n0_vec_list_update[[m]])
        
      }

      
      ### Evaluate loss 
      loss = eval_loss_v4(edge_time_mat_list = edge_time_mat_list,
                          n0_vec_list = n0_vec_list_update,
                          clusters_list = clusters_list,
                          freq_trun = freq_trun,
                          t_vec = t_vec)$loss
      loss_update = loss
    }
    
    if (N_shrink>10) {
      message("[est_n0_vec_v7.1]: Reached maximum shrinkage number.")
    }
    
    
    ### ### ###
    
    ### Save "loss_first_iter"
    if (n_iter==1) {
      loss_first_iter = loss_update
    }    
    
    n0_mat_list_update = n0_vec2mat(n0_vec = n0_vec_list_update)
    
    ### Evaluate stopping criterion
    ### [TO DO: adjust epsilon value]
    converge = (loss_current-loss_update)/(loss_first_iter-loss_update+.Machine$double.eps) < epsilon  
    
    
    ### Debug
    # browser()
    # plot(n0_vec_list_current[[1]], n0_vec_list_update[[1]])
    # step_size*gradient_vec_list[[1]]
    
    
    ### *update -> *current
    n_iter = n_iter + 1
    n0_vec_list_current = n0_vec_list_update
    n0_mat_list_current = n0_mat_list_update
    loss_current = loss_update
    loss_history = c(loss_history, loss_current)
    
  }
  
  if (n_iter>max_iter) {
    message("[est_n0_vec_v7.1]: Reached maximum iteration number.")
  }
  
  ### Debug
  # browser()
  # plot(loss_history[-1])
  # print(n_iter)
  
  
  
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

