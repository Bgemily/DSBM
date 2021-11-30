

# k-means clustering
### Based on v3
### The clustering task considers both conn_prob and shape. No need to specify the weight.
### Use est_n0_vec_v4.1: force n0 to be less than earliest edge time
### Normalize node_cdf_array when updating time shifts.
cluster_kmeans_v6.1 = function(edge_time_mat_list, 
                             clusters_list, n0_vec_list=NULL, n0_mat_list=NULL, center_cdf_array=NULL, 
                             t_vec=seq(0,200,length.out=1000), order_list=NULL, 
                             fix_timeshift=FALSE,
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
                                               n0_mat_list = n0_mat_list, t_vec = t_vec)
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
                                            n0_mat = n0_mat, t_vec = t_vec)
    
    ### Update clusters
    ############### V1: classification. Compare each node with each cluster
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
        dist_array[,q,k] = dist_array_1[,q,k]*(conn_prob_N_vec)^2 + 
                            dist_array_2[,q,k]*(sum((center_cdf_array_normed[q,k,]*
                                                       I(center_cdf_array_normed[q,k,]<1))^2))
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

    ################ V2: clustering. Use pairwise distance, do not need estimated centers
    # pdist_array = array(0, dim=c(N_node, N_node, N_clus))
    # for (l in 1:N_clus) {
    #   pdist_array[ , ,l] = rdist::pdist(node_cdf_array[ ,l, ])
    # }
    # weights = sapply(clusters, length)
    # weights = weights/sum(weights)
    # pdist_mat = apply(pdist_array, c(1,2), function(dij_vec) sqrt(sum(weights*(dij_vec)^2)))
    # membership = cluster::pam(x=pdist_mat, k=N_clus, diss=TRUE, cluster.only=TRUE)
    # clusters = mem2clus(membership)
    # dist_to_centr_vec=NA
    ##############################
    
    clusters_list[[m]] = clusters
  }
  
  
  ### Debug
  # print("####################")
  # print("Forcing the clusters to be true clusters...")
  # print("####################")
  # clusters_list = lapply(1:N_subj,function(.)list(1:30,31:60,61:90))
  #############
  
  ### Update center --- will be the initial value for next step
  center_cdf_array = get_center_cdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                             clusters_list = clusters_list, 
                                             n0_mat_list = n0_mat_list, t_vec = t_vec)
  # browser(); mean(clusters_list[[1]][[1]])
  
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
                                               n0_mat_list = n0_mat_list, t_vec = t_vec)
  } else {
    res = est_n0_vec_v4.1(edge_time_mat_list = edge_time_mat_list, 
                          clusters_list = clusters_list, 
                          n0_vec_list = n0_vec_list, n0_mat_list = n0_mat_list,
                          center_cdf_array = center_cdf_array,
                          t_vec = t_vec, order_list=order_list, ...)
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
# 
# tmp = cluster_kmeans_v2(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list, epsilon=1)
# tmp$clusters_list
# plot(tmp$v_vec_list[[1]], time_shift_list[[1]])


