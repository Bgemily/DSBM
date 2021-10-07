

# transfer vector of membership to list of clusters 
mem2clus = function(membership, N_clus_min=3){
  N_clus = max(membership,N_clus_min)
  # clusters = vector("list", N_clus)
  clusters = list()
  for (l in 1:N_clus) {
    clusters[[l]] = which(membership==l)
  }
  return(clusters)
}


clus2mem = function(clusters){
  membership = unlist(clusters)
  N_clus = length(clusters)
  for (l in 1:N_clus) {
    membership[clusters[[l]]] = l
  }
  return(membership)
}
