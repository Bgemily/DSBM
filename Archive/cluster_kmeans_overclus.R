
# over-clustering (k-means++) and merge (k-medoids)
# Actually, this function doesn't need arguments "clusters" and "n0_vec"

cluster_kmeans_overclus = function(edge_time_mat, clusters, n0_vec, N_overclus, N_clus, MaxIter = 30, N_trial = 10, t_vec=seq(0,50,0.05), bw=1){
  
  # align aggregated pdfs, get aligned_pdf_mat
  N_node = nrow(edge_time_mat)
  node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), n0_vec = numeric(N_node), t_vec = t_vec, bw = bw)
  n0_vec_init = est_n0_vec(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), t_vec = t_vec, bw = bw)
  
  aligned_pdf_mat = t(sapply(1:N_node, function(i)shift(node_pdf_array[i,1,], n0_vec_init[i], pp=TRUE)))
  
  
  # over-cluster shifted_pdf_mat via kmeanspp
  res_overclus = kmeanspp(data_mat = aligned_pdf_mat, N_clus = N_overclus, MaxIter = MaxIter, N_trial = N_trial)
  clusters_overclus = res_overclus$clusters 
  
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters_overclus, t_vec = t_vec, bw = bw)  
  
  # merge clusters
  clusters = merge_clusters(edge_time_mat = edge_time_mat, clusters = clusters_overclus, n0_vec = n0_vec, N_clus = N_clus, t_vec = t_vec, bw = bw)
  
  # udpate n0_vec
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)  
  
  
  # exact-cluster shifted_pdf_mat via kmeanspp
  res_exaclus = kmeanspp(data_mat = aligned_pdf_mat, N_clus = N_clus, MaxIter = MaxIter, N_trial = N_trial)
  clusters_exaclus = res_exaclus$clusters 
  
  
  return(list(clusters=clusters, n0_vec=n0_vec, clusters_exaclus=clusters_exaclus, clusters_overclus=clusters_overclus))
  
}

