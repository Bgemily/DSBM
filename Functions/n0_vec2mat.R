# Turn n0_vec to n0_mat through maximum operations

n0_vec2mat_ = function(n0_vec){
  n0_mat = matrix(nrow=length(n0_vec), ncol=length(n0_vec))
  for (i in 1:nrow(n0_mat)) {
    for (j in 1:ncol(n0_mat)) {
      n0_mat[i,j] = max(n0_vec[i], n0_vec[j])
    }
  }
  return(n0_mat)
}

n0_vec2mat = function(n0_vec){
  if (!is.list(n0_vec)) {
    n0_mat = n0_vec2mat_(n0_vec = n0_vec)
    return(n0_mat)
  }
  else{
    n0_vec_list = n0_vec
    N_subj = length(n0_vec_list)
    n0_mat_list = list()
    for (m in 1:N_subj) {
      n0_mat_list[[m]] = n0_vec2mat_(n0_vec = n0_vec_list[[m]])
    }
    return(n0_mat_list)
  }
  
}


### Test
# n0_vec = 1:10
# n0_vec2mat(n0_vec)
# n0_vec2mat(list(n0_vec-5, n0_vec))
