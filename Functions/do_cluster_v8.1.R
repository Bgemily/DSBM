

### main algorithm
### Use cluster_kmeans_v6.1: Consider both conn_prob and shape(cdf) when clustering.
### Update time shifts at each iteration
do_cluster_v8.1 = function(edge_time_mat_list, N_clus, 
                         clusters_list_init, n0_vec_list_init, n0_mat_list_init,
                         freq_trun=Inf, 
                         total_time = 200, t_vec=seq(0,total_time,length.out=1000),
                         MaxIter=10, conv_thres=5e-3, save_est_history=FALSE, 
                         fix_timeshift=FALSE,
                         ...)
{
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  clusters_history = list()
  loss_history = c()
  cluster_time = align_time = 0
  
  v_vec_history = list()
  center_cdf_array_history = list()

  # Initialize clusters and time shifts -------------------------------------
  
  clusters_list = clusters_list_init
  n0_vec_list = n0_vec_list_init
  n0_mat_list = n0_mat_list_init
  order_list = lapply(n0_vec_list, function(n0_vec)order(n0_vec))
  
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             freq_trun = Inf,
                                             n0_mat_list = n0_mat_list, t_vec = t_vec)
  
  ### Save estimation
  if (save_est_history==TRUE){
    clusters_history = c(clusters_history, list(clusters_list[[1]]))
    v_vec_list = lapply(n0_vec_list, function(vec)vec*t_unit)
    v_vec_history = c(v_vec_history, list(v_vec_list[[1]]))
    center_cdf_array_history= c(center_cdf_array_history, list(center_cdf_array))
  }
  
  
  ### Evaluate loss function
  loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
                      n0_mat_list = n0_mat_list, 
                      clusters_list = clusters_list, 
                      freq_trun = Inf,
                      center_cdf_array = center_cdf_array, 
                      t_vec = t_vec)$loss
  loss_history = c(loss_history, loss)
  
  
  
  # Update clusters and connecting patterns separately for each subject ---------------------

  clusters_list_update = clusters_list_current = clusters_list
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  v_vec_list_update = v_vec_list_current = list()
  center_cdf_array_list = list()
  
  ### Estimate parameters for each subject
  for (m in 1:N_subj) {
    n_iter = 1
    stopping = FALSE
    loss_history_tmp = c()
    
    ### Initialize connecting patterns
    center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list[m], 
                                               clusters_list = clusters_list[m], 
                                               freq_trun = Inf,
                                               n0_mat_list = n0_mat_list[m], t_vec = t_vec)
    center_cdf_array_update = center_cdf_array_current = center_cdf_array
    
    
    while (!stopping & n_iter<=MaxIter){
      ### Update clusters, time shifts and connecting patterns
      res = cluster_kmeans_v6.1(edge_time_mat_list=edge_time_mat_list[m], 
                              clusters_list=clusters_list_current[m], 
                              n0_vec_list=n0_vec_list_current[m], n0_mat_list=n0_mat_list_current[m], 
                              center_cdf_array = center_cdf_array_current,
                              freq_trun = freq_trun,
                              t_vec=t_vec, order_list=NULL, 
                              fix_timeshift=fix_timeshift,
                              ...)
      clusters_list_update[m] = res$clusters_list
      n0_vec_list_update[m] = res$n0_vec_list
      n0_mat_list_update[m] = res$n0_mat_list
      v_vec_list_update[m] = res$v_vec_list
      center_cdf_array_update = res$center_cdf_array
      
      ### Record computing time for clustering and aligning
      cluster_time = cluster_time + res$cluster_time
      align_time = align_time + res$align_time
      
      
      ### Evaluate stopping criterion
      delta_n0_vec = sum((unlist(n0_vec_list_update[m])-unlist(n0_vec_list_current[m]))^2) / 
        ( sum(unlist(n0_vec_list_current[m])^2) + .Machine$double.eps )
      
      
      clusters_update = clusters_list_update[[m]]
      clusters_current = clusters_list_current[[m]]
      delta_clusters = 1 - get_one_ARI(memb_est_vec = clus2mem(clusters_update), 
                                       memb_true_vec = clus2mem(clusters_current))
      
      
      delta_center_cdf = tryCatch(sum((center_cdf_array_update-center_cdf_array_current)^2) / 
                                    (sum(center_cdf_array_current^2) + .Machine$double.eps),
                                  error=function(x)1)
      
      stopping = mean(c(delta_center_cdf,delta_clusters,delta_n0_vec)) < conv_thres
      
      
      
      ### *update -> *current
      n_iter = n_iter+1
      clusters_list_update[m] -> clusters_list_current[m]
      n0_vec_list_update[m] -> n0_vec_list_current[m] 
      n0_mat_list_update[m] -> n0_mat_list_current[m] 
      v_vec_list_update[m] -> v_vec_list_current[m]
      center_cdf_array_update -> center_cdf_array_current 
      
      
      ### Test: Evaluate loss
      loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list[m],
                          n0_mat_list = n0_mat_list_current[m],
                          clusters_list = clusters_list_current[m],
                          center_cdf_array = center_cdf_array_current, 
                          freq_trun = Inf,
                          t_vec = t_vec)$loss
      loss_history = c(loss_history, loss)
      
      ### Save estimation
      if (save_est_history==TRUE) {
        clusters_history = c(clusters_history, list(clusters_list_current[[m]]))
        v_vec_history = c(v_vec_history, list(v_vec_list_current[[m]]))
        center_cdf_array_history= c(center_cdf_array_history, list(center_cdf_array_current))
      }
      
    }
    
    if(n_iter>MaxIter){
      print("[do_cluster_v8.1:] Reached maximum iteration number.")
    }
    N_iteration = n_iter
    
    
    center_cdf_array_list[[m]] = center_cdf_array_current
    
  }
  
  
  # Match clusters across subjects ------------------------------------------
  
  if (N_subj>1) {
    ### Find permutation
    res = match_clusters_v2(center_cdf_array_list = center_cdf_array_list)
    permn_list = res$permn_list
    
    ### Permutate clusters
    for (m in 1:N_subj) {
      permn = permn_list[[m]]
      clusters_list_current[[m]] = clusters_list_current[[m]][permn]
    }
  }
  
  
  center_cdf_array_current = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                     clusters_list = clusters_list_current, 
                                                     n0_mat_list = n0_mat_list_current, 
                                                     freq_trun = Inf,
                                                     t_vec = t_vec)
  
  ### Evaluate loss function
  loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
                      n0_mat_list = n0_mat_list_current, 
                      clusters_list = clusters_list_current, 
                      center_cdf_array = center_cdf_array_current, 
                      freq_trun = Inf,
                      t_vec = t_vec)$loss
  loss_history = c(loss_history, loss)
  
  
  
  # Combine all subjects and update params ----------------------------
  
  if (N_subj>1) {
    n_iter = 1
    stopping = FALSE
    
    ### Update params again
    while (!stopping & n_iter<=MaxIter){
      
      ### Update clusters, time shifts and connecting patterns
      res = cluster_kmeans_v6.1(edge_time_mat_list=edge_time_mat_list, 
                              clusters_list=clusters_list_current, 
                              n0_vec_list=n0_vec_list_current, n0_mat_list=n0_mat_list_current, 
                              center_cdf_array = center_cdf_array_current,
                              t_vec=t_vec, order_list=NULL, 
                              fix_timeshift=fix_timeshift,
                              ...)
      clusters_list_update = res$clusters_list
      n0_vec_list_update = res$n0_vec_list
      n0_mat_list_update = res$n0_mat_list
      v_vec_list_update = res$v_vec_list
      center_cdf_array_update = res$center_cdf_array
      
      ### Record computing time for clustering and aligning
      cluster_time = cluster_time + res$cluster_time
      align_time = align_time + res$align_time
      
      ### Evaluate loss function
      loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
                          n0_mat_list = n0_mat_list_update, 
                          clusters_list = clusters_list_update, 
                          freq_trun = freq_trun,
                          center_cdf_array = center_cdf_array_update, t_vec = t_vec)$loss
      loss_history = c(loss_history, loss)
      
      ### Evaluate stopping criterion
      delta_n0_vec = sum((unlist(n0_vec_list_update)-unlist(n0_vec_list_current))^2) / 
        ( sum(unlist(n0_vec_list_current)^2) + .Machine$double.eps )
      
      delta_clusters_vec = numeric(length = N_subj)
      for (m in 1:N_subj) {
        clusters_update = clusters_list_update[[m]]
        clusters_current = clusters_list_current[[m]]
        delta_tmp = 1 - get_one_ARI(memb_est_vec = clus2mem(clusters_update), 
                                    memb_true_vec = clus2mem(clusters_current))
        delta_clusters_vec[m] = delta_tmp
      }
      weights = N_node_vec/sum(N_node_vec)
      delta_clusters = sum(delta_clusters_vec * weights)
      
      delta_center_cdf = tryCatch(sum((center_cdf_array_update-center_cdf_array_current)^2) / 
                                    (sum(center_cdf_array_current^2) + .Machine$double.eps),
                                  error=function(x)1)
      stopping = mean(c(delta_center_cdf,delta_clusters,delta_n0_vec)) < conv_thres
      
      
      ### *update -> *current
      n_iter = n_iter+1
      clusters_list_update -> clusters_list_current
      n0_vec_list_update -> n0_vec_list_current 
      n0_mat_list_update -> n0_mat_list_current 
      v_vec_list_update -> v_vec_list_current
      center_cdf_array_update -> center_cdf_array_current 
      
    }
    
    
    if (n_iter>MaxIter) {
      message("[do_cluster_v8]: Reached maximum iteration number.")
    }
    
  }
  
  
  
  # Get final result --------------------------------------------------------
  
  clusters_list_current -> clusters_list
  n0_vec_list_current -> n0_vec_list
  n0_mat_list_current -> n0_mat_list
  v_vec_list_current -> v_vec_list
  center_cdf_array_current -> center_cdf_array
  
  ### Save freq_trun_mat 
  if (freq_trun<Inf){
    freq_trun_mat = matrix(data = 0, nrow=N_clus, ncol=N_clus)
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        # N_basis = sum(center_fft_array[q,k,]!=0)
        freq_trun_mat[q,k] = freq_trun
      }
    }
  } else{
    freq_trun_mat = NULL
  }
  
  
  ### Get estimated pdf using kernel smoothing
  center_pdf_array = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, 
                                             t_vec = t_vec)
  
  return(list(clusters_list=clusters_list, 
              v_vec_list=v_vec_list, 
              n0_vec_list=n0_vec_list, n0_mat_list=n0_mat_list,
              center_pdf_array=center_pdf_array, center_cdf_array=center_cdf_array,
              clusters_history=clusters_history, 
              v_vec_history=v_vec_history,
              center_cdf_array_history=center_cdf_array_history,
              loss_history=loss_history,
              freq_trun_mat=freq_trun_mat,
              N_iteration=N_iteration,
              cluster_time=cluster_time, align_time=align_time))
  
}

