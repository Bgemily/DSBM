
### Perform algorithm based on cumulative intensities 
do_cluster_pdf = function(edge_time_mat_list, N_clus, 
                          clusters_list_init, n0_vec_list_init, n0_mat_list_init,
                          freq_trun=15, step_size=200,
                          total_time = 200, t_vec=seq(0,total_time,length.out=1000),
                          MaxIter=10, conv_thres=5e-3, 
                          opt_radius=max(t_vec)/2,
                          fix_timeshift=FALSE,
                          gamma=0.001,
                          ...)
{
  
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  clusters_history = list()
  n0_mat_list_history = list()
  loss_history = c()
  cluster_time = align_time = 0
  
  # Initialize clusters and time shifts -------------------------------------
  
  clusters_list = clusters_list_init
  n0_vec_list = n0_vec_list_init
  n0_mat_list = n0_mat_list_init

  center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                          clusters_list = clusters_list, 
                                          freq_trun = freq_trun, 
                                          n0_mat_list = n0_mat_list, t_vec = t_vec)
  

  ### Save estimation
  clusters_history = c(clusters_history, list(clusters_list))
  n0_mat_list_history = c(n0_mat_list_history, list(n0_mat_list))
  
  
  # Update clusters and connecting patterns separately for each subject ---------------------
  
  
  clusters_list_update = clusters_list_current = clusters_list
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  v_vec_list_update = v_vec_list_current = list()
  center_fft_array_list = list()
  
  ### Estimate parameters for each subject
  for (m in 1:N_subj) {
    n_iter = 1
    stopping = FALSE
    loss_history_tmp = c()
    
    ### Initialize connecting patterns 
    center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list[m], 
                                            clusters_list = clusters_list[m], 
                                            freq_trun = freq_trun, 
                                            n0_mat_list = n0_mat_list[m], t_vec = t_vec)
    center_fft_array_update = center_fft_array_current = center_fft_array
    
    
    while (!stopping & n_iter<=MaxIter){
      
      ### Update clusters, time shifts and connecting patterns 
      res = cluster_kmeans_pdf(edge_time_mat_list=edge_time_mat_list[m], 
                                clusters_list=clusters_list_current[m], 
                                n0_vec_list=n0_vec_list_current[m], n0_mat_list=n0_mat_list_current[m], 
                                center_fft_array = center_fft_array_current,
                                freq_trun = freq_trun, 
                                t_vec=t_vec, 
                                opt_radius=opt_radius,
                                fix_timeshift=fix_timeshift,
                               gamma=gamma,
                                ...)
      clusters_list_update[m] = res$clusters_list
      n0_vec_list_update[m] = res$n0_vec_list
      n0_mat_list_update[m] = res$n0_mat_list
      v_vec_list_update[m] = res$v_vec_list
      center_fft_array_update = res$center_fft_array
      
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
      
      
      delta_center_fft = tryCatch(sum((abs(center_fft_array_update-center_fft_array_current))^2) / 
                                    (sum(abs(center_fft_array_current)^2) + .Machine$double.eps),
                                  error=function(x)1)
      
      stopping = mean(c(delta_center_fft,delta_clusters,delta_n0_vec)) < conv_thres
      
      
      
      ### *update -> *current
      n_iter = n_iter+1
      clusters_list_update[m] -> clusters_list_current[m]
      n0_vec_list_update[m] -> n0_vec_list_current[m] 
      n0_mat_list_update[m] -> n0_mat_list_current[m] 
      v_vec_list_update[m] -> v_vec_list_current[m]
      center_fft_array_update -> center_fft_array_current 
      
      
      clusters_history = c(clusters_history, list(clusters_list_current[m]))
      n0_mat_list_history = c(n0_mat_list_history, list(n0_mat_list_current[m]))
      
      
    }
    
    if (n_iter>MaxIter) {
      message("[do_cluster_pdf]: Reached maximum iteration number.")
    }
    N_iteration = n_iter
    
    
    center_fft_array_list[[m]] = center_fft_array_current
    
  }
  
  
  # Match clusters across subjects ------------------------------------------
  
  if (N_subj>1) {
    ### Find permutation [TO DO: finish to-do's in find_permn_v3.]
    res = match_clusters_v4(center_fft_array_list = center_fft_array_list)
    permn_list = res$permn_list
    
    ### Permutate clusters
    for (m in 1:N_subj) {
      permn = permn_list[[m]]
      clusters_list_current[[m]] = clusters_list_current[[m]][permn]
    }
  }
  
  
  
  
  
  clusters_history = c(clusters_history, list(clusters_list_current))
  
  
  center_fft_array_current = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                                  clusters_list = clusters_list_current, 
                                                  n0_mat_list = n0_mat_list_current, 
                                                  freq_trun = freq_trun,  t_vec = t_vec)
  
  clusters_history = c(clusters_history, list(clusters_list_current))
  n0_mat_list_history = c(n0_mat_list_history, list(n0_mat_list_current))
  
  
  
  # Combine all subjects and update params ----------------------------
  
  if (N_subj>1) {
    n_iter = 1
    stopping = FALSE
    
    ### Update params again
    while (!stopping & n_iter<=MaxIter){
      
      ### Update clusters, time shifts and connecting patterns
      res = cluster_kmeans_pdf(edge_time_mat_list=edge_time_mat_list, 
                                clusters_list=clusters_list_current, 
                                n0_vec_list=n0_vec_list_current, n0_mat_list=n0_mat_list_current, 
                                center_pdf_array = center_pdf_array_current,
                                t_vec=t_vec, 
                                opt_radius=opt_radius,
                                fix_timeshift=fix_timeshift,
                               gamma=gamma,
                                ...)
      clusters_list_update = res$clusters_list
      n0_vec_list_update = res$n0_vec_list
      n0_mat_list_update = res$n0_mat_list
      v_vec_list_update = res$v_vec_list
      center_pdf_array_update = res$center_pdf_array
      
      ### Record computing time for clustering and aligning
      cluster_time = cluster_time + res$cluster_time
      align_time = align_time + res$align_time
      
      
      clusters_history = c(clusters_history, list(clusters_list_update))
      

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
      
      delta_center_pdf = tryCatch(sum((center_pdf_array_update-center_pdf_array_current)^2) / 
                                    (sum(center_pdf_array_current^2) + .Machine$double.eps),
                                  error=function(x)1)
      stopping = mean(c(delta_center_pdf,delta_clusters,delta_n0_vec)) < conv_thres
      
      
      ### *update -> *current
      n_iter = n_iter+1
      clusters_list_update -> clusters_list_current
      n0_vec_list_update -> n0_vec_list_current 
      n0_mat_list_update -> n0_mat_list_current 
      v_vec_list_update -> v_vec_list_current
      center_pdf_array_update -> center_pdf_array_current 
      
    }
    
    
    if (n_iter>MaxIter) {
      message("[do_cluster_pdf]: Reached maximum iteration number.")
    }
    
  }
  
  
  
  # Get final result --------------------------------------------------------
  
  clusters_list_current -> clusters_list
  n0_vec_list_current -> n0_vec_list
  n0_mat_list_current -> n0_mat_list
  v_vec_list_current -> v_vec_list
  center_fft_array_current -> center_fft_array
  

  center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                          clusters_list = clusters_list, 
                                          freq_trun = freq_trun, 
                                          n0_mat_list = n0_mat_list, t_vec = t_vec)
  ### Get freq_trun_mat from center_fft_array
  freq_trun_mat = matrix(data = 0, nrow=N_clus, ncol=N_clus)
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      # N_basis = sum(center_fft_array[q,k,]!=0)
      freq_trun_mat[q,k] = freq_trun
    }
  }
  
  ### Convert fourier series back to the (smoothed) pdf
  center_pdf_array = array(dim=c(N_clus, N_clus, length(t_vec)))
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      fft_truncated = center_fft_array[q,k,]
      func = Re(fft(c(tail(fft_truncated, freq_trun+1), 
                      rep(0, length(t_vec)-2*freq_trun-1),
                      head(fft_truncated, freq_trun)), inverse = TRUE))
      center_pdf_array[q,k,] = func
    }
  }
  

  return(list(clusters_list=clusters_list, 
              loss_history=loss_history,
              clusters_history=clusters_history, n0_mat_list_history=n0_mat_list_history,
              v_vec_list=v_vec_list,
              center_pdf_array=center_pdf_array, 
              center_fft_array=center_pdf_array,
              freq_trun_mat=freq_trun_mat,
              N_iteration=N_iteration,
              cluster_time=cluster_time, align_time=align_time))
  
}


