# update n0_vec
est_n0_vec = function(edge_time_mat, clusters, center_cdf_array=NULL, L_mat=NULL, 
                      t_vec=seq(0,50,0.05), bw=1, standardize=FALSE, step_size=0.02){
  
  # step_size = 0.02
  
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(edge_time_mat)
  
  n0_vec = numeric(N_node)
  
  # do not shift any event at begining
  ########### V1
  # node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = numeric(N_node), n0_mat = matrix(0,N_node,N_node),
                                      # t_vec = t_vec, bw = bw)
  ############ V2
  # node_cdf_array = get_node_cdf_array(edge_time_mat = edge_time_mat, clusters = clusters, 
  #                                     n0_mat = 0, t_vec = t_vec, standardize=standardize)
  ############ V3
  earliest_edge_time = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
  order = order(earliest_edge_time) # order patterning time, use this order as order of v_i's
  edge_time_mat_tmp = edge_time_mat
  for (i in 1:(length(order)-1)) {
    edge_time_mat_tmp[order[i],order[(i+1):length(order)]] = NA #only consider {j:v_j<=v_i}
  }
  node_cdf_array = get_node_cdf_array(edge_time_mat = edge_time_mat_tmp, clusters = clusters, 
                                      n0_mat = 0, t_vec = t_vec, standardize=standardize)
  #############################
  
  
  
  for (q in 1:length(clusters)) {
    if (length(clusters[[q]]) == 1) {
      n0_vec[clusters[[q]]] = 0
      next
    }
    
    ########## V1
    # node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat, clusters = clusters, n0_vec = numeric(N_node), n0_mat = matrix(0,N_node,N_node),
    #                                     t_vec = t_vec, bw = bw)
    # # align each node with the first node (i_0) in this cluster
    # i_0 = clusters[[q]][1]
    # for (i in clusters[[q]]) {
    #   # align i-th row towards i_0-th row
    #   possible_n0 = get_dist_betw_pdfarray(pdf_array_1 = node_pdf_array[i, , , drop=F],
    #                                        pdf_array_2 = node_pdf_array[i_0, , , drop=F],
    #                                        symmetric=FALSE, t_unit = t_unit, t_vec = t_vec)$n0_mat
    #   if(length(clusters)>1) possible_n0 = possible_n0[-q] # only consider time shifts in connection with other clusters. This is crucial.
    #   n0_vec[i] = possible_n0[which.max(abs(possible_n0))]
    # }
    
    ########### V2
    if (!is.null(center_cdf_array)){
      # align each node in this cluster with the center of this cluster
      for (i in clusters[[q]]) {
        f1_list = lapply(1:length(clusters), function(l) node_cdf_array[i,l, ])
        f2_list = lapply(1:length(clusters), function(l) center_cdf_array[q,l, ])
        
        ##### V1
        f1_list = lapply(1:length(clusters), function(l) f1_list[[l]]*
                           max(f2_list[[l]])/max(max(f1_list[[l]]), 1e-6))
        ####################
        
        ##### V1
        # weights = L_mat[q, ] * sapply(clusters, function(clus) length(setdiff(clus, i))) # L_{q,l}*|Gamma_l\{i}|
        ##### V2
        # weights = sapply(clusters, function(clus) length(intersect(clus, order[1:i]))) # |{j: v_j<=v_i && z_j==l}|
        ##### V3
        weights = sapply(clusters, function(clus) length(intersect(clus, order[1:which(order==i)]))) # |{j: v_j<=v_i && z_j==l}|
        ##################
        if(sum(weights)==0)
          weights = NULL
        else
          weights = weights/sum(weights)
        
        n0_vec[i] = align_multi_curves_gd(f1_list = f1_list, f2_list = f2_list,
                                          step_size = step_size, pp = FALSE, 
                                          t_unit = t_unit, weights = weights)$n0
      }
    }
    else{
      print("center_cdf_array is NULL.")
      # align each node with the first node (i_0) in this cluster
      i_0 = clusters[[q]][1]
      for (i in clusters[[q]]) {
        # align i-th row towards i_0-th row
        f1_list = lapply(1:length(clusters), function(l) node_cdf_array[i,l, ])
        f2_list = lapply(1:length(clusters), function(l) node_cdf_array[i_0,l, ])
        weights = sapply(clusters, function(clus) length(setdiff(clus, i))) # |Gamma_l\{i}|
        if(sum(weights)==0)
          weights = NULL
        else
          weights = weights/sum(weights)

        n0_vec[i] = align_multi_curves_gd(f1_list = f1_list, f2_list = f2_list,
                                          step_size = step_size, pp = FALSE, t_unit = t_unit, weights = weights)$n0
      }
    }
    
    ##################################
    
    
    ############ V1
    # make the minimum time shift to be zero
    # n0_vec[clusters[[q]]] = n0_vec[clusters[[q]]] - min(n0_vec[clusters[[q]]])
    ############ V2
    # n0_vec[clusters[[q]]] = n0_vec[clusters[[q]]] - min(n0_vec[clusters[[q]]])
    # earliest_edge_time = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
    # n0_vec[clusters[[q]]] = n0_vec[clusters[[q]]] + min(earliest_edge_time[clusters[[q]]])/t_unit
    ############ V3
    # n0_vec[clusters[[q]]] = (earliest_edge_time[clusters[[q]]])/t_unit
    ############ V4
    # make the mean time shift to be zero
    # n0_vec[clusters[[q]]] = n0_vec[clusters[[q]]] - mean(n0_vec[clusters[[q]]])
    ############ V5
    # let min == earliest edge time
    earliest_edge_time = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
    cluswise_min = sapply(seq(length(clusters)), function(k)min(earliest_edge_time[clusters[[k]]]))
    clus_id = which(cluswise_min <= max(earliest_edge_time[clusters[[q]]]))
    min_n0 = min(edge_time_mat[clusters[[q]], unlist(clusters[clus_id])]
                 [edge_time_mat[clusters[[q]], unlist(clusters[clus_id])]>1])/t_unit
    n0_vec[clusters[[q]]] = n0_vec[clusters[[q]]] - min(n0_vec[clusters[[q]]]) + min_n0
    #############################
    
  }
  
  ##### V1
  n0_vec = n0_vec - min(n0_vec)
  ##### V3
  # n0_vec = n0_vec - 10/t_unit
  ##### V2
  # earliest_edge_time = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
  # if(sum(n0_vec>(earliest_edge_time/t_unit))>0){
  #   n0_vec = n0_vec - max(n0_vec-earliest_edge_time/t_unit)
  # }
  ####################

  return(n0_vec)
}

