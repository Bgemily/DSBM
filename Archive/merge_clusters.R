# obtain center_pdf_array and merge similar clusters

merge_clusters = function(edge_time_mat, clusters, n0_vec, n0_mat=NULL, N_clus, t_vec=seq(0,50,0.05), bw=1, intensity=TRUE){
  N_overclus = length(clusters)
  if(N_overclus <= N_clus)
    return(clusters)
  
  # get pairwise distance between clusters
  center_pdf_array = get_center_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
                                          t_vec = t_vec, bw = bw)
  clus_degree_mat = get_clus_degree_mat(edge_time_mat = edge_time_mat, clusters = clusters, intensity=intensity)
  dist_mat = pairwise_dist_mat(pdf_array = center_pdf_array, degree_mat = clus_degree_mat, t_vec = t_vec)$dist_mat
  
  # cluster over-specified clusters via k-medoids
  clus_membership = cluster::pam(x = dist_mat, k = N_clus)$cluster
  
  # merge clusters
  node_membership_overclus = clus2mem(clusters)
  membership = numeric(nrow(edge_time_mat))
  for (i in 1:length(membership)) {
    membership[i] = clus_membership[node_membership_overclus[i]]
  }
  clusters = mem2clus(membership)
  
  return(clusters)
}