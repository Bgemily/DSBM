# Obtain aggregated pdf_vec list (N_node*1)
get_aggre_pdf_list = function(edge_time_mat, pairwise_dist, dist_thres, t_vec, tau_mat, true_pdf_fun_list, membership_true){
  aggre_pdf_list = list() 
  N_node = nrow(edge_time_mat)
  weight_mat = I(pairwise_dist<dist_thres) / rowSums(pairwise_dist<dist_thres)
  for (i in 1:N_node) {
    tmp_mat = sapply(which(pairwise_dist[i,]<dist_thres), function(j) sapply(t_vec-tau_mat[i,j], true_pdf_fun_list[[membership_true[i]]][[membership_true[j]]]))
    aggre_pdf_list[[i]] = rowMeans(tmp_mat)
  }
  return(aggre_pdf_list)
}
