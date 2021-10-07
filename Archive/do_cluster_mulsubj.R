

# main algorithm (with combining multiple subjects)
do_cluster_mulsubj = function(edge_time_mat_list, N_clus, N_overclus=N_clus, MaxIter=3, t_vec=seq(0, 50, 0.05), bw=1){
  
  N_node = nrow(edge_time_mat_list[[1]])
  
  # initial clustering (spectral clustering)
  clusters = list(c(1:N_node))
  n0_mat = matrix(0, N_node, N_node)
  
  res_list = lapply(edge_time_mat_list, function(x) cluster_specc_overclus(x, clusters, n0_mat, N_overclus, N_clus, t_vec, bw))
  
  clusters_list = lapply(res_list, '[[', 'clusters')
  n0_mat_list = lapply(res_list, '[[', 'n0_mat')
  
  # save the clustering result of spectral clustering
  clusters_list_specc = clusters_list
  n0_mat_list_specc = n0_mat_list
  clusters_specc_exaclus = lapply(res_list, '[[', 'clusters_exaclus')

    
  # k-means
  for (. in 1:MaxIter) {
    n0_vec_list = mapply(get_n0_vec, n0_mat_list, clusters_list, SIMPLIFY = FALSE) # n0_vec_list: invariant with permutations of clusters
    
    # match clusters so that *clusters are of the same order* across subjects
    clusters_list = match_clusters(edge_time_mat_list=edge_time_mat_list, clusters_list=clusters_list, 
                                   n0_vec_list=n0_vec_list, t_vec=t_vec, bw=bw)$clusters_list
    
    # obtain connecting patterns between clusters by combining subjects
    center_pdf_array = get_center_pdf_array_mulsubj(edge_time_mat_list, clusters_list, n0_vec_list, t_vec, bw)
    
    res_list = mapply(cluster_kmeans, edge_time_mat_list, clusters_list, n0_mat_list, center_pdf_array=list(center_pdf_array), t_vec=list(t_vec), bw=list(bw), SIMPLIFY = FALSE)    
    
    clusters_list = lapply(res_list, '[[', 'clusters')
    n0_mat_list = lapply(res_list, '[[', 'n0_mat')
    
  }
  
  
  # get center_pdf_array_list
  n0_vec_list = mapply(get_n0_vec, n0_mat_list, clusters_list, SIMPLIFY = FALSE)
  center_pdf_array_list = mapply(get_center_pdf_array, edge_time_mat_list, clusters_list, n0_vec_list, SIMPLIFY = FALSE)

    
  
  # return lists of clusters and n0_mat
  return(list(clusters=clusters_list, n0_mat=n0_mat_list, center_pdf_array_list=center_pdf_array_list, clusters_specc=clusters_list_specc, n0_mat_specc=n0_mat_list_specc, clusters_specc_exaclus=clusters_specc_exaclus))
  
}

