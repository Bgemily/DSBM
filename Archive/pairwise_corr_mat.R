

pairwise_corr_mat = function(pdf_array, degree_mat, t_vec=seq(0, 50, 0.05)){
  res = pairwise_dist_mat(pdf_array = pdf_array, degree_mat = degree_mat, t_unit = NULL, t_vec = t_vec)
  dist_mat = res$dist_mat
  corr_mat = exp(-dist_mat^2/median(dist_mat)^2)
  return(list(corr_mat = corr_mat, n0_mat = res$n0_mat, dist_mat = dist_mat))
}

