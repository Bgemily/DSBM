

# main algorithm
do_cluster = function(edge_time_mat, N_clus, N_overclus=N_clus, MaxIter=10, N_trial=10, 
                      t_vec=seq(0, 200, length.out=1000), bw=1, standardize=FALSE, step_size=0.02, conv_thres=1e-3){
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(edge_time_mat)
  clusters_history = list()
  cluster_time = align_time = 0
  
  # initialize clusters and n0_mat
  res = get_init(edge_time_mat=edge_time_mat, N_clus=N_overclus, 
                 MaxIter=MaxIter, N_trial=N_trial, t_vec=t_vec, bw=bw, 
                 standardize=standardize)
  clusters = res$clusters
  n0_mat = res$n0_mat
  n0_vec = res$n0_vec
  L_mat = res$L_mat

  clusters_history = c(clusters_history, list(clusters))
  
  center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_mat = n0_mat, t_vec = t_vec)
  
  
  # record merged clusters
  if(N_overclus>N_clus){
    clusters_merged = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                     N_clus = N_clus, t_vec = t_vec, bw = bw)
    clusters_history = c(clusters_history, list(clusters_merged))
  }
  
  
  if(N_overclus>N_clus){
    # continue over-clustering
    for (. in 1:3){
      ############ V1
      # res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec, n0_mat=n0_mat,
      #                      t_vec=t_vec, bw=bw, center_pdf_array = NULL)
      # clusters = res$clusters
      # n0_vec = res$n0_vec
      # n0_mat = res$n0_mat
      ############ V2
      res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_mat=n0_mat, 
                           L_mat = L_mat, t_vec=t_vec, step_size = step_size)
      clusters = res$clusters
      n0_vec = res$n0_vec
      n0_mat = res$n0_mat
      L_mat = res$L_mat
      ##############################
      
      clusters_history = c(clusters_history, list(clusters))
      
      # merge clusters
      if (N_overclus>N_clus){
        clusters_merged = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                         N_clus = N_clus, t_vec = t_vec, bw = bw)
        clusters_history = c(clusters_history, list(clusters_merged))
      }
      
    }
  }
  
  
  
  # merge clusters
  if (N_overclus>N_clus){
    ############### V1
    clusters = clusters_merged
    res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
    n0_mat = res$n0_mat
    n0_vec = res$n0_vec
    ############### V2
    # res = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters,...)
    # clusters = res$clusters
    # n0_mat = res$n0_mat
    # n0_vec = res$n0_vec
    # L_mat = res$L_mat
    ###################################
  }
  
  

  # continue k-means
  clusters_old = NULL
  n_iter = 1
  
  stopping = FALSE
  Loss_c = Inf
  loss_history = c()
  while (!stopping & n_iter<=MaxIter){
    clusters_old = clusters
    n0_vec_old = n0_vec
    center_cdf_array_old = center_cdf_array
    ############## V1
    # res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec, n0_mat=n0_mat,
    #                      t_vec=t_vec, bw=bw, center_pdf_array = NULL)
    # clusters = res$clusters
    # n0_vec = res$n0_vec
    # n0_mat = res$n0_mat
    ############# V2
    res = cluster_kmeans(edge_time_mat=edge_time_mat, clusters=clusters, n0_mat=n0_mat, 
                         L_mat = L_mat, t_vec=t_vec, standardize=standardize, step_size = step_size)
    clusters = res$clusters
    n0_vec = res$n0_vec
    n0_mat = res$n0_mat
    L_mat = res$L_mat
    cluster_time = cluster_time + res$cluster_time
    align_time = align_time + res$align_time
    ################################

    clusters_history = c(clusters_history, list(clusters))
    
    ########### V1
    # stopping = identical(clusters, clusters_old)
    # loss_history = NA
    ########### V2
    # Loss_u = eval_loss(edge_time_mat = edge_time_mat, n0_mat = n0_mat, clusters = clusters, 
    #                    t_vec = t_vec, standardize = FALSE)$loss
    # stopping = abs(Loss_c-Loss_u)/Loss_u < conv_thres
    # Loss_c = Loss_u
    # loss_history = c(loss_history, Loss_u)
    ############ V3
    center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_mat = n0_mat, t_vec = t_vec)
    
    delta_n0_vec = sum((n0_vec_old-n0_vec)^2) / sum(n0_vec_old^2)
    delta_clusters = 1 - get_one_ARI(clus2mem(clusters_old), clus2mem(clusters))
    delta_center_cdf = tryCatch(sum((center_cdf_array_old-center_cdf_array)^2) / sum(center_cdf_array^2),
                                error=function(x)1)
      
    
    stopping = mean(c(delta_center_cdf,delta_clusters,delta_n0_vec)) < conv_thres
    ##########################
    
    n_iter = n_iter+1
  }


  
  # # or stick with spectral clustering
  # clusters_history = list(clusters)
  # for (. in 1:MaxIter){
  #   res = cluster_specc(edge_time_mat, clusters, n0_vec, N_clus, t_vec, bw)
  #   clusters = res$clusters
  #   n0_vec = res$n0_vec
  #   clusters_history = c(clusters_history, list(clusters))
  # }
  
  
  # get center_pdf_array
  center_pdf_array = get_center_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                          n0_vec = n0_vec, n0_mat = n0_mat, t_vec = t_vec, bw = bw)
  center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_mat = n0_mat, t_vec = t_vec)
  
  
  return(list(clusters=clusters, clusters_history=clusters_history, n0_vec=n0_vec, n0_mat=n0_mat, L_mat=L_mat, 
              center_pdf_array=center_pdf_array, center_cdf_array=center_cdf_array,
              cluster_time=cluster_time, align_time=align_time, loss_history=loss_history))
  
}

