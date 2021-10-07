


get_ICL = function(edge_time_mat_list, clusters_list, v_vec_list, center_pdf_array){
  
  N_subj = length(edge_time_mat_list)
  N_clus = length(clusters_list[[1]])
  N_node_vec = sapply(edge_time_mat_list, nrow)
  N_node = sum(N_node_vec)

  max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
  t_vec = seq(0, max_time+10, length.out = 1000)
  v_mat_list = n0_vec2mat(v_vec_list)
  n0_mat_list = lapply(v_mat_list, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
  
  center_cdf_array = center_pdf_array
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      center_cdf_array[q,k,] = cumsum(center_pdf_array[q,k,])
    }
  }
  
  loss = 0
  for (m in 1:N_subj) {
    node_pdf_array = get_node_pdf_array(edge_time_mat = edge_time_mat_list[[m]], 
                                           clusters = clusters_list[[m]], 
                                           n0_mat = n0_mat_list[[m]], t_vec = t_vec, bw=10)
    
    for (q in 1:N_clus) {
      for (i in clusters_list[[m]][[q]]) {
        for (k in 1:N_clus) {
          node_pdf_vec = node_pdf_array[i,k,]
          center_pdf_vec = center_pdf_array[q,k,]
          
          loss_tmp = sum((node_pdf_vec-center_pdf_vec)^2)
          N_edge = sum(edge_time_mat_list[[m]][i,clusters_list[[m]][[k]]]<Inf)
          
          if (!is.na(loss_tmp*N_edge)) {
            loss = loss + loss_tmp*N_edge
          }
          else{
            loss = loss
          }
          
        }
      }
    }
  }
  
  
  ICL = -log(loss) - 1/2*(N_clus-1)*log(N_node) - 1/2*( log(sum(N_node_vec*(N_node_vec-1)))*
                                                          (N_clus*(N_clus+1)/2) )
  ICL = -log(loss)
  return(ICL)
}




# Choose N_clus for real data ---------------------------------------------------------------


# ICL_vec = c()
# file_list = list.files(path = '../Results/Rdata/real_data_results/', 
#                        pattern = 'Partial_subj_N_clus*', full.names = TRUE)
# for (file in file_list) {
#   load(file)
#   ICL = get_ICL(edge_time_mat_list = edge_time_mat_list[avai_subj], 
#                 clusters_list = clusters_list, 
#                 v_vec_list = v_vec_list,
#                 center_pdf_array = center_pdf_array)
#   ICL_vec = c(ICL_vec, ICL)
# }
# plot(2:5,ICL_vec, type='b')
