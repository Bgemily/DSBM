# update n0_mat

est_n0_mat_v2 = function(edge_time_mat_list, 
                         clusters_list, center_cdf_array=NULL, 
                         n0_vec_list=NULL, n0_mat_list=NULL,
                         t_vec=seq(0,50,0.05), step_size=0.02)
{
  
  n0_vec_list = est_n0_vec_v2(edge_time_mat_list = edge_time_mat_list, 
                              clusters_list = clusters_list, center_cdf_array=center_cdf_array, 
                              n0_vec_list = n0_vec_list, n0_mat_list = n0_mat_list,
                              t_vec = t_vec, step_size=step_size)
  
  
  n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  
  return(list(n0_mat_list=n0_mat_list, n0_vec_list=n0_vec_list))
  
}

