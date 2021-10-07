

# k-means clustering
cluster_kmeans = function(edge_time_mat, clusters, n0_vec=NULL, n0_mat=NULL, L_mat=NULL, 
                          center_cdf_array=NULL, t_vec=seq(0, 50, 0.05), bw=1, intensity=TRUE, 
                          standardize=FALSE, step_size=0.02){
  
  t_unit = t_vec[2]-t_vec[1]
  N_node = nrow(edge_time_mat); 
  N_clus = length(clusters)
  
  if(!is.null(center_cdf_array) && length(clusters)!=dim(center_cdf_array)[2]) 
    stop("dim(center_pdf_array)[2] and length(clusters) should match.")
  else if(!is.null(center_cdf_array) && length(clusters)!=dim(center_cdf_array)[1])
    N_clus = dim(center_cdf_array)[1]
  
  
  # update clusters--------------------------------------------------------------------------
  
  cluster_start_time = Sys.time()
  ############### V1
  # # estimate node_pdf_array
  # node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
  #                                     t_vec = t_vec, bw = bw)
  # 
  # # update center_pdf_array
  # if (is.null(center_pdf_array)) {
  #   center_pdf_array = get_center_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = n0_vec, n0_mat = n0_mat,
  #                                           t_vec = t_vec, bw = bw)
  # }
  
  ################ V2
  # estimate node_cdf_array
  node_cdf_array = get_node_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                      n0_mat = n0_mat, t_vec = t_vec, standardize=standardize)
  
  # update center_cdf_array
  if (is.null(center_cdf_array)) {
    center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
                                            n0_mat = n0_mat, t_vec = t_vec, standardize=standardize)
  }
  ############################################
  
  
  # update clusters
  ############### V1
  # membership = numeric(N_node)
  # dist_to_centr_vec = numeric(N_node)
  # # degree_mat = get_node_degree_mat(edge_time_mat = edge_time_mat, clusters = clusters, intensity=intensity)
  # for (i in 1:N_node) {
  #   weights = sapply(clusters, function(clus) length(setdiff(clus,i)))
  #   weights = weights/sum(weights)
  # 
  #   dist_vec = numeric(N_clus)
  #   for (l in 1:N_clus) {
  #     # distance between ith node and l-th cluster
  #     dist_vec[l] = get_dist_betw_pdfarray(pdf_array_1 = node_cdf_array[i, , , drop=F],
  #                                          pdf_array_2 = center_cdf_array[l, , ,drop=F],
  #                                          symmetric=FALSE, weights=weights, t_unit=t_unit, t_vec = t_vec,
  #                                          use_shift_inv_dist = F, pp=FALSE)$dist
  #   }
  #   membership[i] = which.min(dist_vec)
  #   dist_to_centr_vec[i] = min(dist_vec)
  # }
  # 
  # clusters = mem2clus(membership)
  ################ V2: use pairwise distance, do not need estimated centers
  pdist_array = array(0, dim=c(N_node, N_node, N_clus))
  for (l in 1:N_clus) {
    pdist_array[ , ,l] = rdist::pdist(node_cdf_array[ ,l, ])
  }
  weights = sapply(clusters, length)
  weights = weights/sum(weights)
  pdist_mat = apply(pdist_array, c(1,2), function(dij_vec) sqrt(sum(weights*(dij_vec)^2)))
  res = cluster::pam(x=pdist_mat, k=N_clus, diss=TRUE, cluster.only=TRUE)
  clusters = mem2clus(res)
  # update center
  center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters,
                                          n0_mat = n0_mat, t_vec = t_vec, standardize=standardize)
  dist_to_centr_vec=NA
  ##############################
  
  cluster_end_time = Sys.time()
  cluster_time = cluster_end_time - cluster_start_time
  
  
  # update n0_mat------------------------------------------------------------------------------
  
  align_start_time = Sys.time()
  ########## V1
  # # n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  # res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  # n0_mat = res$n0_mat
  # n0_vec = res$n0_vec
  ########## V2
  res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, center_cdf_array = center_cdf_array,
                   L_mat = L_mat, t_vec = t_vec, standardize=standardize, step_size=step_size)
  n0_mat = res$n0_mat
  n0_vec = res$n0_vec
  L_mat = res$L_mat
  ########## V3
  # center_cdf_array = get_center_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters,
  #                                         n0_mat = n0_mat, t_vec = t_vec, standardize=TRUE)
  # print("Use standardization when estimating n0_vec.")
  # res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, center_cdf_array = center_cdf_array,
  #                  L_mat = L_mat, t_vec = t_vec, standardize=TRUE)
  # n0_mat = res$n0_mat
  # n0_vec = res$n0_vec
  # L_mat = res$L_mat
  ###########################
  align_end_time = Sys.time()
  align_time = align_end_time - align_start_time
  
  
  return(list(clusters=clusters, n0_vec=n0_vec, n0_mat=n0_mat, L_mat=L_mat, dist_to_centr_vec=dist_to_centr_vec,
              cluster_time=cluster_time, align_time=align_time))
}

