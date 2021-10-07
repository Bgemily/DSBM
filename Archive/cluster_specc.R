

# clustering
cluster_specc = function(edge_time_mat, clusters, n0_vec, N_clus, t_vec=seq(0,50,0.05), bw=1){
  
  pdf_array = get_node_pdf_array( edge_time_mat=edge_time_mat, clusters=clusters, n0_vec=n0_vec, t_vec=t_vec, bw=bw)
  degree_mat = get_node_degree_mat(edge_time_mat, clusters, intensity=TRUE)
  
  res = pairwise_corr_mat(pdf_array, degree_mat = degree_mat)
  
  membership = spectral_clustering(res$corr_mat, N_clus)
  clusters = mem2clus(membership)
  
  # n0_mat = res$n0_mat
  dist_mat = res$dist_mat
  corr_mat = res$corr_mat
  
  n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  
  return(list(clusters=clusters, n0_vec=n0_vec, dist_mat=dist_mat, corr_mat=corr_mat))
}

