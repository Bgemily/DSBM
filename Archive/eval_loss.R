### Loss = Sum of |S^{-v}N-F|^2.

eval_loss = function(edge_time_mat_list, n0_vec_list, clusters_list, 
                     t_vec=seq(0,50,0.05)){

  N_subj = length(edge_time_mat_list)
  N_clus = length(clusters_list[[1]])
  t_unit = t_vec[2]-t_vec[1]
  
  
  ### Evaluate loss (up to a constant)
  loss = 0
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      adjus_node_cdf_bigmat = c()
      n0_mat_longvec = c()
      for (m in 1:N_subj) {
        ### V1
        adjus_edge_time_mat = edge_time_mat_list[[m]]-n0_vec2mat(n0_vec_list[[m]])*t_unit
        ### v2
        # adjus_edge_time_mat = edge_time_mat_list[[m]]-n0_vec2mat(n0_vec_list[[m]])*t_unit*0
        ### ### ###
        clusters = clusters_list[[m]]
        adjus_node_cdf_mat = t(sapply(adjus_edge_time_mat[clusters[[q]], clusters[[k]]],
                                      function(t) if(!is.na(t)) ecdf(t)(t_vec)
                                      else rep(NA,length(t_vec))) )
        non_NA_id = which(!is.na(rowSums(adjus_node_cdf_mat)))
        adjus_node_cdf_mat = adjus_node_cdf_mat[non_NA_id, ]
        adjus_node_cdf_bigmat = rbind(adjus_node_cdf_bigmat, adjus_node_cdf_mat)
        
        ### n0_mat_longvec: length \approx M*n^2
        n0_mat_vec = c(n0_vec2mat(n0_vec_list[[m]]))
        n0_mat_longvec = c(n0_mat_longvec, n0_mat_vec)
        n0_mat_longvec = round(n0_mat_longvec)
      }
      ### Sweep out col_means
      loss_qk = sweep(adjus_node_cdf_bigmat, 2,
                      apply(adjus_node_cdf_bigmat, 2, mean))
      ### Shift back towards right
      for (i in 1:nrow(loss_qk)) {
        loss_qk[i,] = c( numeric(n0_mat_longvec[i]), loss_qk[i, 1:(ncol(loss_qk)-n0_mat_longvec[i])])
      } 
      
      
      loss_qk = sum(loss_qk^2)
      loss = loss + loss_qk
    }
  }
  
  return(list(loss=loss))
}





# res = results1$conn_prob_0.2[[1]]
# eval_loss(edge_time_mat = res$network$edge_time_mat,n0_mat = res$clus_result$n0_mat,
#           clusters=res$clus_result$clusters,t_vec = res$network$t_vec,standardize = F)
# eval_loss(edge_time_mat = res$network$edge_time_mat,n0_mat = res$network$tau_mat/res$network$t_vec[2],
#           clusters=mem2clus(res$network$membership_true),t_vec = res$network$t_vec,standardize = F)

