

### main algorithm
### For real data analysis
### Based on v11
### est_n0_vec_v5 ==> est_n0_vec_v5.1, match_clusters ==> match_clusters_v3
### For stage I&II: cluster_kmeans_v6.1
 
do_cluster_v11.1 = function(edge_time_mat_list, N_clus, 
                         clusters_list_init, n0_vec_list_init, n0_mat_list_init,
                         total_time = 200, t_vec=seq(0,total_time,length.out=1000),
                         MaxIter=10, conv_thres=5e-3, bw=5, ...)
{
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  clusters_history = list()
  loss_history = c()
  cluster_time = align_time = 0
  
# Initialize clusters and time shifts -------------------------------------
  
  clusters_list = clusters_list_init
  n0_vec_list = n0_vec_list_init
  n0_mat_list = n0_mat_list_init
  order_list = lapply(n0_vec_list, function(n0_vec)order(n0_vec))
 
  clusters_history = c(clusters_history, list(clusters_list))
  
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, t_vec = t_vec)
  
  
  
  ### Evaluate loss function
  loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
                      n0_mat_list = n0_mat_list, 
                      clusters_list = clusters_list, 
                      center_cdf_array = center_cdf_array, t_vec = t_vec)$loss
  loss_history = c(loss_history, loss)
  
 

# Update clusters and connecting patterns separately for each subject ---------------------

  
  clusters_list_update = clusters_list_current = clusters_list
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  n0_mat_list_update = n0_mat_list_current = n0_mat_list
  v_vec_list_update = v_vec_list_current = list()
  center_cdf_array_list = list()
  center_pdf_array_list = list()
  
  ### Estimate parameters for each subject
  for (m in 1:N_subj) {
    n_iter = 1
    stopping = FALSE
    loss_history_tmp = c()
    
    ### Initialize connecting patterns
    center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list[m], 
                                               clusters_list = clusters_list[m], 
                                               n0_mat_list = n0_mat_list[m], t_vec = t_vec)
    center_cdf_array_update = center_cdf_array_current = center_cdf_array
    
    
    while (!stopping & n_iter<=MaxIter){
      
      ### Update clusters, time shifts and connecting patterns
      res = cluster_kmeans_v6.1(edge_time_mat_list=edge_time_mat_list[m],
                              clusters_list=clusters_list_current[m],
                              n0_vec_list=n0_vec_list_current[m], n0_mat_list=n0_mat_list_current[m],
                              center_cdf_array = center_cdf_array_current,
                              t_vec=t_vec, order_list=NULL, ...)
      
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
      # loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list[m],
      #                     n0_mat_list = n0_mat_list_current[m],
      #                     clusters_list = clusters_list_current[m],
      #                     center_cdf_array = center_cdf_array_current, t_vec = t_vec)$loss
      # loss_history_tmp = c(loss_history_tmp, loss)
      
    }
    
  
    center_cdf_array_list[[m]] = center_cdf_array_current
    
    ### Compute center_pdf_array for each individual
    center_pdf_array = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list[m], 
                                               clusters_list = clusters_list_current[m], 
                                               n0_mat_list = n0_mat_list_current[m],
                                               t_vec = t_vec, bw=bw)
    center_pdf_array_list[[m]] = center_pdf_array
    
  }
  
  order_list_current = lapply(n0_vec_list_current, function(n0_vec)order(n0_vec))
  
  ### Test
  # clusters_list_current -> clusters_list
  # n0_vec_list_current -> n0_vec_list
  # n0_mat_list_current -> n0_mat_list
  # v_vec_list_current -> v_vec_list
  # 
  # return(list(clusters_list=clusters_list, clusters_history=clusters_history, 
  #             v_vec_list=v_vec_list,
  #             center_pdf_array_list = center_pdf_array_list))
  
  

# Match clusters across subjects ------------------------------------------

  if (N_subj>1) {
    ### Find permutation
    res = match_clusters_v2(center_cdf_array_list = center_cdf_array_list)
    permn_list = res$permn_list
    
    ### Permute clusters
    for (m in 1:N_subj) {
      permn = permn_list[[m]]
      clusters_list_current[[m]] = clusters_list_current[[m]][permn]
    }
    
    ### Permute center_pdf_array_list 
    for (m in 1:N_subj) {
      permn = permn_list[[m]]
      center_pdf_array_list[[m]] = center_pdf_array_list[[m]][permn, permn, ]
    }
    
  }
  
  
  
  clusters_history = c(clusters_history, list(clusters_list_current))
  
  
  center_cdf_array_current = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                     clusters_list = clusters_list_current, 
                                                     n0_mat_list = n0_mat_list_current, t_vec = t_vec)
  
  
  ### Evaluate loss function
  loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
                      n0_mat_list = n0_mat_list_current, 
                      clusters_list = clusters_list_current, 
                      center_cdf_array = center_cdf_array_current, t_vec = t_vec)$loss
  loss_history = c(loss_history, loss)
  
  

# Combine all subjects and update params ----------------------------

  ### Update time shifts based on merged connection patterns
  res = est_n0_vec_v5.1(edge_time_mat_list = edge_time_mat_list, 
                      clusters_list = clusters_list_current, 
                      n0_vec_list = n0_vec_list_current, n0_mat_list = n0_mat_list_current,
                      center_cdf_array = center_cdf_array_current,
                      t_vec = t_vec, order_list=NULL, ...)
  n0_vec_list_current = res$n0_vec_list
  n0_mat_list_current = res$n0_mat_list
  center_cdf_array_current = res$center_cdf_array
  
  
  ### Update clusters, conn_patt, and conn_prob alternatively
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
                              t_vec=t_vec, order_list=NULL, ...)
      clusters_list_update = res$clusters_list
      n0_vec_list_update = res$n0_vec_list
      n0_mat_list_update = res$n0_mat_list
      v_vec_list_update = res$v_vec_list
      center_cdf_array_update = res$center_cdf_array
      
      ### Record computing time for clustering and aligning
      cluster_time = cluster_time + res$cluster_time
      align_time = align_time + res$align_time
      
      
      clusters_history = c(clusters_history, list(clusters_list_update))
      
      
      ### Evaluate loss function
      # loss = eval_loss_v2(edge_time_mat_list = edge_time_mat_list, 
      #                     n0_mat_list = n0_mat_list_update, 
      #                     clusters_list = clusters_list_update, 
      #                     center_cdf_array = center_cdf_array_update, t_vec = t_vec)$loss
      # loss_history = c(loss_history, loss)
      
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
      message("[do_cluster_v11]: Reached maximum iteration number.")
    }
    
  }
  
  

# Get final result --------------------------------------------------------

  clusters_list_current -> clusters_list
  n0_vec_list_current -> n0_vec_list
  n0_mat_list_current -> n0_mat_list
  v_vec_list_current -> v_vec_list
  center_cdf_array_current -> center_cdf_array
  
  center_pdf_array = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, t_vec = t_vec, bw=bw)
  
  for (m in 1:N_subj) {
    center_pdf_array_list[[m]] = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list[m], 
                                                         clusters_list = clusters_list[m], 
                                                         n0_mat_list = n0_mat_list[m], 
                                                         t_vec = t_vec, bw=bw)
  }
  
  return(list(clusters_list=clusters_list, clusters_history=clusters_history, 
              loss_history=loss_history,
              v_vec_list=v_vec_list,
              center_pdf_array=center_pdf_array, center_cdf_array=center_cdf_array,
              center_pdf_array_list = center_pdf_array_list,
              cluster_time=cluster_time, align_time=align_time))
  
}


# Test --------------------------------------------------------------------

# res = generate_network2_v2(SEED=831, N_subj=3, N_node_vec=rep(90,3),
#                            N_clus=3,
#                            total_time=200,
#                            conn_patt_var=1, conn_patt_sep = 1.3, const=40, conn_prob_mean = 1, conn_prob_rad = 0,
#                            time_shift_struc=max, time_shift_rad = 10)
# 
# edge_time_mat_list = res$edge_time_mat_list
# cdf_true_array = res$cdf_true_array
# pdf_true_array = res$pdf_true_array
# clusters_list = res$clus_true_list
# time_shift_list = res$time_shift_list
# 
# tmp = do_cluster_v2(edge_time_mat_list = edge_time_mat_list, N_clus = 3, max_iter=1, step_size=0.02)
# plot(tmp$loss_history)
# tmp$clusters_list
# tmp$clusters_history
# plot(time_shift_list[[1]][], tmp$v_vec_list[[1]][])

# tmp2 = do_cluster(edge_time_mat = edge_time_mat_list[[1]], N_clus = 3)
# tmp2$clusters
# tmp2$clusters_history[[1]]
# plot(time_shift_list[[1]], tmp2$n0_vec)
