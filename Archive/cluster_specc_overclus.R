
# over-clustering and merge (spectral clustering)
cluster_specc_overclus = function(edge_time_mat, clusters, n0_mat, N_overclus, N_clus, t_vec=seq(0,50,0.05), bw=1){
  
  if (N_overclus==N_clus) return(cluster_specc(edge_time_mat, clusters, n0_mat, N_clus, t_vec, bw))
  
  
  # over-clustering
  res = cluster_specc(edge_time_mat, clusters, n0_mat, N_overclus, t_vec, bw)
  clusters_overclus = res$clusters
  n0_mat = res$n0_mat
  
  # save results of exact-clustering
  corr_mat = res$corr_mat
  membership_exaclus = spectral_clustering(corr_mat, N_clus)
  clusters_exaclus = mem2clus(membership_exaclus)

  
  # obtain centers' pdf_array
  n0_vec = get_n0_vec(n0_mat, clusters_overclus)
  center_pdf_array = get_center_pdf_array(edge_time_mat, clusters_overclus, n0_vec, t_vec, bw)
  
  clus_degree_mat = get_clus_degree_mat(edge_time_mat, clusters_overclus, intensity=TRUE)
  
  # merge clusters
  res = pairwise_corr_mat(center_pdf_array, degree_mat = clus_degree_mat)
  clus_membership = spectral_clustering(res$corr_mat, N_clus) # res$corr_mat: N_overclus * N_overclus
  
  node_membership_overclus = clus2mem(clusters_overclus)
  
  membership = numeric(nrow(edge_time_mat))
  for (i in 1:length(membership)) {
    membership[i] = clus_membership[node_membership_overclus[i]]
  }
  clusters = mem2clus(membership)
  
  return(list(clusters=clusters, n0_mat=n0_mat, clusters_exaclus=clusters_exaclus, corr_mat=corr_mat))
}

