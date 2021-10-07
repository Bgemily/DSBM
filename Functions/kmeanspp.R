# implementation of k-means++

kmeanspp = function(data_mat, N_clus, MaxIter = 100, N_trial = 10){
  N_node = nrow(data_mat)
  
  out = list()
  out$tot.withinss = Inf
  
  
  for (trial in 1:N_trial){
    center_ids = numeric(N_clus)
    center_ids[1] = sample(N_node, 1)
    
    for (l in 2:N_clus) {
      sq_dist_to_centr_mat = apply(data_mat[center_ids, , drop=F], 1, function(center)rowSums((sweep(data_mat, 2, center))^2))
      sq_dist_to_centr_vec = apply(sq_dist_to_centr_mat, 1, min)
      sq_dist_to_centr_vec[center_ids] = 0
      center_ids[l] = sample(N_node, 1, prob = sq_dist_to_centr_vec)
    }
    
    tmp_out = kmeans(data_mat, centers = data_mat[center_ids, , drop=F], iter.max = MaxIter)
    
    tmp_out$init_center_id = center_ids
    
    if (tmp_out$tot.withinss < out$tot.withinss) {
      out <- tmp_out
    }
  }
  
  out$membership = out$cluster
  out$cluster = NULL
  out$clusters = mem2clus(out$membership)
  
  return(out)
}

# test
# init_kmeans(shifted_pdf_mat, 5)$clusters


