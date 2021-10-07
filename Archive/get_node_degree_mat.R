
# get degrees from each node to each cluster
get_node_degree_mat = function(edge_time_mat, clusters, intensity=TRUE){
  degree_mat = matrix(0, nrow=nrow(edge_time_mat), ncol=length(clusters))
  # 
  # if(intensity){
  #   for (i in 1:nrow(degree_mat)) {
  #     for (l in 1:ncol(degree_mat)) {
  #       degree_mat[i,l] = length(clusters[[l]])
  #     }
  #   } 
  #   return(degree_mat)
  # }
  # 
  # else{
    for (i in 1:nrow(degree_mat)) {
      for (l in 1:ncol(degree_mat)) {
        degree_mat[i,l] = sum(edge_time_mat[i, clusters[[l]]]<Inf)
      }
    } 
    return(degree_mat)
  # }
  
}

