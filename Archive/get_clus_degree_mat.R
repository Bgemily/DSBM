
# get degrees from each cluster to another cluster
get_clus_degree_mat = function(edge_time_mat, clusters, intensity=TRUE){
  degree_mat = matrix(0, nrow=length(clusters), ncol=length(clusters))
  
  # if(intensity){
  #   for (q in 1:nrow(degree_mat)) {
  #     for (l in 1:ncol(degree_mat)) {
  #       degree_mat[q,l] = length(clusters[[q]]) * length(clustes[[l]])
  #     }
  #   } 
  #   return(degree_mat)
  # }
  # 
  # else{
    for (q in 1:nrow(degree_mat)) {
      for (l in 1:ncol(degree_mat)) {
        degree_mat[q,l] = sum(edge_time_mat[clusters[[q]], clusters[[l]]]<Inf)
        if(q==l) degree_mat[q,l] = degree_mat[q,l] %/% 2 # because each edge is counted twice
      }
    } 
    return(degree_mat)
  # }
  
}

