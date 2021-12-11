
# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# -------------------------------------------------------------------------
load("../Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v11_freqtrun7/pr=0.9,n=30,beta=1.3,V=80/beta/1.3/N_trial10_20211206_231811.Rdata")

for (res in results) {
  network_param = res$network_param
  # Generate networks 
  
  network_list = do.call(what = generate_network2_v3, args = network_param)
  
  edge_time_mat_list = network_list$edge_time_mat_list
  pdf_true_array = network_list$pdf_true_array
  membership_true_list = network_list$membership_true_list
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$time_shift_list
  order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))
  
  center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list,
                                                 clusters_list = res$clusters_list_est,
                                                 n0_mat_list = list(n0_vec2mat(res$v_vec_list_est)[[1]]/(network_param$t_vec[2])),
                                                 t_vec = network_param$t_vec)
  print(max(center_pdf_array_est[,,1]))
  
}

print(max(center_pdf_array_est))



# -------------------------------------------------------------------------


