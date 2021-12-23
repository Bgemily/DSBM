
# pdf_array:  n*k*N, degree_mat: n*k
pairwise_dist_mat = function(pdf_array, degree_mat, t_unit = 0.05, t_vec = seq(0, 50, 0.05)){
  N_node = dim(pdf_array)[1]
  dist_mat = matrix(0, N_node, N_node)
  n0_mat = matrix(0, N_node, N_node)
  for (i in 1:(N_node-1)) {
    for (j in ((i+1):N_node)) {
      # NEED JUSTIFICATION
      weights = degree_mat[i,]*degree_mat[j,]; weights = sqrt(weights); weights = weights/sum(weights)
      
      # weights = NULL
      
      res = get_dist_betw_pdfarray(pdf_array_1 = pdf_array[i, , , drop=F], pdf_array_2 = pdf_array[j, , , drop=F], 
                                   symmetric=FALSE, weights=weights, t_unit = t_unit, t_vec = t_vec)

      dist_mat[i,j] = res$dist
      dist_mat[j,i] = dist_mat[i,j]
      
      n0 = res$n0_mat[which.max(abs(res$n0_mat))] ### NEED JUSTIFICATION
      n0_mat[i,j] = n0
      n0_mat[j,i] = -n0_mat[i,j]
    }
  }
  return(list(dist_mat=dist_mat, n0_mat=n0_mat))
}

