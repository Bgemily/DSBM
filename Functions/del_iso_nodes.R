
#  delete isolated nodes
del_iso_nodes = function(network){
  edge_time_mat = network$edge_time_mat
  
  i=1
  while (i<=nrow(edge_time_mat)) {
    if(sum(edge_time_mat[i,]<Inf)==0)
    {
      cat('Deleting node', i, '\n')
      edge_time_mat = edge_time_mat[-i,-i]
      network$edge_time_mat = network$edge_time_mat[-i, -i]
      network$nodes_mat = network$nodes_mat[-i, ]
      network$tau_vec = network$tau_vec[-i]
      network$pdf_list = network$pdf_list[-i]
      network$membership_true = network$membership_true[-i]
    }
    else
      i = i+1
  }
  return(network)
}

