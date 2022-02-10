

# get shifts of nodes from n0_mat
get_n0_vec = function(n0_mat, clusters){
  n0_vec = vector('numeric', length = nrow(n0_mat))
  for (l in 1:length(clusters)) {
    n0_vec[clusters[[l]]] = n0_mat[ clusters[[l]], clusters[[l]][1] ]
    n0_vec[clusters[[l]]] = n0_vec[clusters[[l]]] - min(n0_vec[clusters[[l]]])
  }
  return(n0_vec)
}

