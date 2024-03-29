

# Perform centering step once and clustering step once
cluster_kmeans_cdf = function(edge_time_mat_list, 
                             clusters_list, n0_vec_list=NULL, n0_mat_list=NULL, 
                             center_cdf_array=NULL, 
                             freq_trun = Inf,
                             t_vec=seq(0,200,length.out=1000), 
                             fix_timeshift=FALSE,
                             gamma=0.1,
                             ...)
{
  
  t_unit = t_vec[2]-t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  
  if (length(unique(sapply(clusters_list,length)))>1) {
    stop("The number of clusters is not the same across subjects.")
  }
  N_clus = length(clusters_list[[1]])
  
  if(!is.null(center_cdf_array) && N_clus!=dim(center_cdf_array)[2]) 
    stop("dim(center_pdf_array)[2] and N_clus do not match.")
  
  if (is.null(n0_vec_list)) {
    res = get_init_v2(edge_time_mat_list=edge_time_mat_list, N_clus=N_clus, t_vec=t_vec)
    n0_vec_list = res$n0_vec_list
    n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  }
  
  
  # Update clusters--------------------------------------------------------------------------
  
  clustering_start_time = Sys.time()
  
  ### Compute center_cdf_array if not given
  if (is.null(center_cdf_array)) {
    center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                               clusters_list = clusters_list, 
                                               n0_mat_list = n0_mat_list, 
                                               freq_trun = Inf, 
                                               t_vec = t_vec)
  }
  
  
  ### Update clusters for all subjects
  for (m in 1:N_subj) {
    N_node = N_node_vec[m]
    edge_time_mat = edge_time_mat_list[[m]]
    clusters = clusters_list[[m]]
    n0_mat = n0_mat_list[[m]]
    
    ### Compute node_cdf_array
    node_cdf_array = get_node_cdf_array_v2(edge_time_mat = edge_time_mat, 
                                            clusters = clusters, 
                                            n0_mat = n0_mat, 
                                           freq_trun = Inf, 
                                           t_vec = t_vec)
    
    ### Update clusters: Compare each node with each cluster
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
    node_cdf_array_normed = node_cdf_array
    center_cdf_array_normed = center_cdf_array
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        ### Compute distance between connecting probabilities
        conn_prob_N_vec = node_cdf_array[,k,dim(node_cdf_array)[3]]
        conn_prob_F = max(center_cdf_array[q,k,])
        dist_array_2[,q,k] = (conn_prob_N_vec - conn_prob_F)^2
        
        ### Normalize (match connecting probabilities)
        node_cdf_array_normed[,k,] = node_cdf_array[,k,] / (node_cdf_array[,k,dim(node_cdf_array)[3]] + 1e-10)
        center_cdf_array_normed[q,k,] = center_cdf_array[q,k,] / max(max(center_cdf_array[q,k,]), 1e-6)
        
        ### Compute distance between normalized distribution
        tmp = rowSums( ( node_cdf_array_normed[,k,] - matrix(data=center_cdf_array_normed[q,k,],
                                                    nrow=N_node, ncol=length(t_vec), byrow = TRUE) )^2 )
        dist_array_1[,q,k] = tmp
        
        ### Weighted sum of conn_prob and shape
        dist_array[,q,k] = dist_array_1[,q,k]*(conn_prob_N_vec) + dist_array_2[,q,k]*gamma*length(t_vec)
        
      }
    }
    
    ### Compute distances
    for (q in 1:N_clus) {
      dist_mat[,q] = rowSums(dist_array[,q,]*weights)
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
  
  ### Update center --- will be the initial value for next step
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, 
                                             freq_trun = Inf, 
                                             t_vec = t_vec)
  
  clustering_end_time = Sys.time()
  cluster_time = clustering_end_time - clustering_start_time
  
  
  # Update time shifts ------------------------------------------------------------------------------
  
  align_start_time = Sys.time()
  if (fix_timeshift) {
    n0_vec_list = n0_vec_list
    n0_mat_list = n0_mat_list
    v_vec_list = lapply(n0_vec_list, function(n0_vec)n0_vec*t_unit)
    v_mat_list = lapply(n0_mat_list, function(n0_mat)n0_mat*t_unit)
    center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                               clusters_list = clusters_list, 
                                               n0_mat_list = n0_mat_list, 
                                               freq_trun = Inf, 
                                               t_vec = t_vec)
  } else {
    res = est_n0_vec_cdf(edge_time_mat_list = edge_time_mat_list, 
                          clusters_list = clusters_list, 
                          n0_vec_list = n0_vec_list, n0_mat_list = n0_mat_list,
                          center_cdf_array = center_cdf_array,
                          freq_trun = Inf, 
                          t_vec = t_vec, ...)
    n0_vec_list = res$n0_vec_list
    n0_mat_list = res$n0_mat_list
    v_vec_list = res$v_vec_list
    v_mat_list = res$v_mat_list
    center_cdf_array = res$center_cdf_array
  }
  
  
  align_end_time = Sys.time()
  align_time = align_end_time - align_start_time

  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list=clusters_list, 
              n0_vec_list=n0_vec_list, n0_mat_list=n0_mat_list, 
              v_vec_list=v_vec_list, v_mat_list=v_mat_list,
              center_cdf_array=center_cdf_array,
              cluster_time=cluster_time, align_time=align_time))
}


