

# Perform centering step once and clustering step once
cluster_kmeans_pdf = function(edge_time_mat_list, clusters_list, 
                             n0_vec_list=NULL, n0_mat_list=NULL, 
                             center_fft_array=NULL, 
                             freq_trun=5, 
                             t_vec=seq(0,200,length.out=1000), 
                             opt_radius=max(t_vec)/2,
                             fix_timeshift=FALSE,
                             gamma=0.001,
                             ...)
{
  
  t_unit = t_vec[2]-t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  if (length(unique(sapply(clusters_list,length)))>1) {
    stop("The number of clusters is not the same across subjects.")
  }
  N_clus = length(clusters_list[[1]])
  
  if(!is.null(center_fft_array) && N_clus!=dim(center_fft_array)[2]) 
    stop("dim(center_pdf_array)[2] and N_clus do not match.")
  
  if (is.null(n0_vec_list)) {
    res = get_init_v2(edge_time_mat_list=edge_time_mat_list, N_clus=N_clus, t_vec=t_vec)
    n0_vec_list = res$n0_vec_list
    n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  }
  
  
  # Update clusters--------------------------------------------------------------------------
  
  clustering_start_time = Sys.time()
  
  ### Compute center_fft_array if not given
  if (is.null(center_fft_array)) {
    center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                            clusters_list = clusters_list, 
                                            freq_trun = freq_trun, 
                                            n0_mat_list = n0_mat_list, 
                                            t_vec = t_vec)
  }
  center_fft_array_normed = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                                 clusters_list = clusters_list, 
                                                 freq_trun = freq_trun, 
                                                 n0_mat_list = n0_mat_list, 
                                                 rmv_conn_prob = TRUE,
                                                 t_vec = t_vec)
  
  
  
  ### Update clusters for all subjects
  for (m in 1:N_subj) {
    N_node = N_node_vec[m]
    edge_time_mat = edge_time_mat_list[[m]]
    clusters = clusters_list[[m]]
    n0_mat = n0_mat_list[[m]]
    
    ### Compute node_fft_array: N_node * N_clus * (2*freq_trun+1)
    node_fft_array = get_node_fft_array(edge_time_mat = edge_time_mat, 
                                        clusters = clusters, 
                                        n0_mat = n0_mat, 
                                        t_vec = t_vec, 
                                        freq_trun = freq_trun)
    node_fft_array_normed = get_node_fft_array(edge_time_mat = edge_time_mat, 
                                               clusters = clusters, 
                                               n0_mat = n0_mat, 
                                               t_vec = t_vec, 
                                               rmv_conn_prob = TRUE, 
                                               freq_trun = freq_trun)
    
    ### Update clusters. Compare each node with each cluster
    membership = numeric(N_node)
    dist_to_centr_vec = numeric(N_node)

    ### Compute weights due to varying cluster sizes
    weights = matrix(nrow=N_node, ncol=N_clus)
    flag_obs_mat = matrix(data=1, nrow=N_node, ncol=N_node) - diag(nrow=N_node, ncol=N_node)
    for (q in 1:N_clus) {
      N_obs = rowSums(as.matrix(flag_obs_mat[,clusters[[q]]]))
      weights[,q] = N_obs / (N_node-1)
    }

    ### Compute distance between each node and each cluster
    dist_mat = matrix(nrow=N_node, ncol=N_clus)
    dist_array = dist_array_1 = dist_array_2 = array(dim=c(N_node,N_clus,N_clus))
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        
        ### Compute distance between all nodes and cluster q w.r.t. clus k
        dist_array_1[,q,k] = length(t_vec) * rowSums( ( abs( node_fft_array_normed[,k,] - 
                                                           matrix(data=center_fft_array_normed[q,k,],
                                                                  nrow=N_node, 
                                                                  ncol=length(center_fft_array_normed[q,k,]), 
                                                                  byrow = TRUE) ) )^2 )
        conn_prob_N_vec = rowSums(edge_time_mat[,clusters[[k]],drop=FALSE]<Inf) / 
          ( rep(length(clusters[[k]]), nrow(edge_time_mat)) - 1*(1:nrow(edge_time_mat) %in% clusters[[k]]) )
        conn_prob_F = sum(edge_time_mat[clusters[[q]], clusters[[k]]]<Inf) / 
          ( length(edge_time_mat[clusters[[q]], clusters[[k]]]) - I(q==k)*length(clusters[[q]]) )
        dist_array_2[,q,k] = (conn_prob_N_vec - conn_prob_F)^2
        
        dist_array[,q,k] = dist_array_1[,q,k]*(conn_prob_N_vec) + dist_array_2[,q,k]*gamma*length(t_vec)
        
      }
    }
    
    ### Compute distances between each node and each clus
    for (q in 1:N_clus) {
      dist_mat[,q] = rowSums(dist_array[,q,]*weights,na.rm = TRUE)
    }
    
      

    ### Update memberships and clusters
    for (i in 1:N_node) {
      dist_vec = dist_mat[i, ]
      membership[i] = which.min(dist_vec)
      dist_to_centr_vec[i] = min(dist_vec)
    }
    clusters = mem2clus(membership = membership, N_clus_min = N_clus)

    clusters_list[[m]] = clusters
  }
  
  
  clustering_end_time = Sys.time()
  cluster_time = clustering_end_time - clustering_start_time
  
  
  # Update time shifts ------------------------------------------------------------------------------
  
  align_start_time = Sys.time()
  
  if (fix_timeshift) {
    n0_vec_list = n0_vec_list
    n0_mat_list = n0_mat_list
    v_vec_list = lapply(n0_vec_list, function(n0_vec)n0_vec*t_unit)
    v_mat_list = lapply(n0_mat_list, function(n0_mat)n0_mat*t_unit)
    center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                            clusters_list = clusters_list, 
                                            freq_trun = freq_trun, 
                                            n0_mat_list = n0_mat_list, 
                                            t_vec = t_vec)
  } else {
    res = est_n0_vec_pdf(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list,
                            n0_vec_list = n0_vec_list, n0_mat_list = n0_mat_list, 
                            freq_trun = freq_trun, 
                            t_vec = t_vec, 
                            opt_radius=opt_radius,
                            ...)
    n0_vec_list = res$n0_vec_list
    n0_mat_list = res$n0_mat_list
    v_vec_list = res$v_vec_list
    v_mat_list = res$v_mat_list
    center_fft_array = res$center_fft_array
    
  }
  

  align_end_time = Sys.time()
  align_time = align_end_time - align_start_time

  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list=clusters_list, 
              n0_vec_list=n0_vec_list, n0_mat_list=n0_mat_list, 
              v_vec_list=v_vec_list, v_mat_list=v_mat_list,
              center_fft_array=center_fft_array,
              cluster_time=cluster_time, align_time=align_time))
}


