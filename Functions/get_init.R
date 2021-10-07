# initialization of clusters and n0_mat
get_init = function(edge_time_mat, N_clus, MaxIter=10, N_trial=10, t_vec=seq(0, 50, 0.05), bw=1, standardize=FALSE){
  time_unit = t_vec[2] - t_vec[1]
  N_node = nrow(edge_time_mat)
  
  # initialize clusters via kmeanspp, initialize n0_mat 
  node_cdf_array = get_node_cdf_array(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), 
                                      n0_mat = matrix(0,N_node,N_node), t_vec = t_vec, standardize=standardize)
  ########### V1
  # n0_vec = est_n0_vec(edge_time_mat = edge_time_mat, clusters = list(c(1:N_node)), t_vec = t_vec, bw = bw)
  ########### V2
  earliest_edge_time = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
  n0_vec = (earliest_edge_time)/time_unit
  n0_vec = round(n0_vec)
  ############ V3
  # earliest_edge_time = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
  # n0_vec = (earliest_edge_time)/time_unit*0
  # n0_vec = round(n0_vec)
  ###################################
  
  aligned_cdf_mat = t(sapply(1:N_node, function(i)shift(node_cdf_array[i,1,], n0_vec[i])))
  
  ######## V1
  # res = kmeanspp(data_mat = aligned_cdf_mat, N_clus = N_clus, MaxIter = MaxIter, N_trial = N_trial)
  # clusters = res$clusters
  ######### V2
  res = cluster::pam(x=aligned_cdf_mat, k=N_clus, diss=FALSE, cluster.only=TRUE)
  clusters = mem2clus(res)
  ############
  
  ####### V1
  # res = est_n0_mat(edge_time_mat = edge_time_mat, clusters = clusters, t_vec = t_vec, bw = bw)
  # n0_mat = res$n0_mat
  # n0_vec = res$n0_vec
  ###########
  
  ###### V2
  # for (l in 1:length(clusters)) {
  #   n0_vec[clusters[[l]]] = n0_vec[clusters[[l]]] - min(n0_vec[clusters[[l]]])
  # }
  # 
  # # initialize L_mat and n0_mat
  # L_mat = matrix(data = NA, nrow = N_clus, ncol = N_clus)
  # n0_mat = matrix(data = NA, nrow = N_node, ncol = N_node)
  # for (q in 1:N_clus) {
  #   for (l in 1:N_clus) {
  #     edge_time_submat = edge_time_mat[clusters[[q]], clusters[[l]], drop=F]
  # 
  # 
  #     tau_q_vec = time_unit * n0_vec[clusters[[q]]]
  #     tau_l_vec = time_unit * n0_vec[clusters[[l]]]
  # 
  # 
  #     edge_time_submat_1 = sweep(edge_time_submat, 1, tau_q_vec) # shift according to cluster q
  #     var_1 = var(edge_time_submat_1[is.finite(edge_time_submat_1)])
  # 
  #     edge_time_submat_2 = sweep(edge_time_submat, 2, tau_l_vec) # shift according to cluster l
  #     var_2 = var(edge_time_submat_2[is.finite(edge_time_submat_2)])
  # 
  # 
  #     if(var_1<var_2||is.na(var_1)) {
  #       # edge_time_submat = edge_time_submat_1
  #       n0_mat[clusters[[q]], clusters[[l]]] = matrix(n0_vec[clusters[[q]]], nrow = length(clusters[[q]]), ncol = length(clusters[[l]]))
  #       n0_mat[clusters[[l]], clusters[[q]]] = t(n0_mat[clusters[[q]], clusters[[l]]])
  #       L_mat[q,l] = 1; L_mat[l,q] = 0; if(q==l) L_mat[q,l] = 1
  #     }
  #     else {
  #       # edge_time_submat = edge_time_submat_2
  #       n0_mat[clusters[[q]], clusters[[l]]] = matrix(n0_vec[clusters[[l]]], nrow = length(clusters[[q]]), ncol = length(clusters[[l]]), byrow = TRUE)
  #       n0_mat[clusters[[l]], clusters[[q]]] = t(n0_mat[clusters[[q]], clusters[[l]]])
  #       L_mat[q,l] = 0; L_mat[l,q] = 1; if(q==l) L_mat[q,l] = 1
  #     }
  # 
  #     if(q==l){
  #       for (i in clusters[[q]]) {
  #         for (j in clusters[[l]]) {
  #           n0_mat[i,j] = min(n0_vec[i], n0_vec[j])
  #         }
  #       }
  #     }
  # 
  #   }
  # }
  # 
  ########## V3
  n0_vec = n0_vec - min(n0_vec)
  L_mat = matrix(data = 1, nrow = N_clus, ncol = N_clus)
  n0_mat = matrix(data = NA, nrow = N_node, ncol = N_node)
  for (i in 1:N_node) {
    for (j in 1:N_node) {
      n0_mat[i,j] = max(n0_vec[i], n0_vec[j])
    }
  }
  ##################
  
  return(list(clusters=clusters, n0_mat=n0_mat, n0_vec=n0_vec, L_mat=L_mat))
}