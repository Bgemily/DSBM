# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Load simulation result and get network parameters -----------------------

load("../Results/Rdata/SNR_Vnot0_v3/main_v5_cdf_v1/pr=0.9,n=30,beta=1.3,V=80/beta/1.9/N_trial10_20211123_211833.Rdata")
network_param = results[[2]]$network_param

# Generate networks -------------------------------------------------------

network_list = do.call(what = generate_network2_v3, args = network_param)

edge_time_mat_list = network_list$edge_time_mat_list
pdf_true_array = network_list$pdf_true_array
membership_true_list = network_list$membership_true_list
clus_true_list = network_list$clus_true_list
v_true_list = network_list$time_shift_list
order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))


# Set algorithm parameters ------------------------------------------------
N_clus = network_param$N_clus
N_clus_min=N_clus 
N_clus_max=N_clus

freq_trun=15 
bw=5 
conv_thres=1e-2 
MaxIter=5
jitter_time_rad = 10 
max_iter=10
total_time = network_param$total_time
opt_radius=total_time/2

t_vec = network_param$t_vec
N_subj = network_param$N_subj


# Apply our method --------------------------------------------------------

res_list = list()
for (N_clus_tmp in N_clus_min:N_clus_max) {
  ### Get initialization -----------
  res = get_init_v4(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp, 
                    t_vec = t_vec)
  
  clusters_list_init = res$clusters_list
  n0_vec_list_init = res$n0_vec_list
  n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
  
  ### Evaluate accuracy of initial time shifts
  spearman_corr_vec = numeric(length = N_subj)
  for (m in 1:N_subj) {
    spearman_corr_vec[m] = cor(v_true_list[[m]], n0_vec_list_init[[m]], method = "spearman")
  }
  spearman_corr_vec_mean_init = mean(spearman_corr_vec)
  
  
  
  # Apply algorithm ---------
  
  
  ### Estimation z,v,f based on cdf
  res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                        total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                        clusters_list_init = clusters_list_init,
                        n0_vec_list_init = n0_vec_list_init, n0_mat_list_init = n0_mat_list_init,
                        )
  
  # Save results of N_clus_tmp ----------------------------------------------
  
  res_list = c(res_list, list(res))
  
  # ### Given true time shifts, update estimations
  # clusters_list_init = res$clusters_list
  # n0_vec_list_init = lapply(v_true_list, function(vec) round(vec/(t_vec[2]-t_vec[1])))
  # n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
  # res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
  #                       total_time = total_time, max_iter=max_iter, t_vec=t_vec,
  #                       clusters_list_init = clusters_list_init,
  #                       n0_vec_list_init = n0_vec_list_init, n0_mat_list_init = n0_mat_list_init)
  # 
  # # Save results of N_clus_tmp ----------------------------------------------
  # 
  # res_list = c(res_list, list(res))
  # 
}

# Plot estimated intensities ----------------------------------------------
plot(y=results[[1]]$v_vec_list[[1]], x=v_true_list[[1]]); abline(a=0,b=1,col='red')

permn = c(1,3,2)
plot_pdf_array_v2(pdf_array_list = list(res_list[[1]]$center_pdf_array[permn,permn,]), 
                  pdf_true_array = pdf_true_array, 
                  clus_size_vec = sapply(res_list[[1]]$clusters_list[[1]][permn],length),
                  y_lim = c(0,0.04),
                  t_vec = t_vec)

