
# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# -------------------------------------------------------------------------

network_list = do.call(generate_network, results[[9]]$network_param)

edge_time_mat_list = network_list$edge_time_mat_list # Generated connecting time matrix


freq_trun = Inf # Cut-off frequency (corresponding to \ell_0 in the main text). For cumulative-intensity-based algorithm, set freq_trun to be Inf.
N_clus_tmp = 3 # Desired number of clusters 
MaxIter = 10 # Maximal number of iterations between the centering and aligning steps
conv_thres = 0.01 # Threshold in convergence criterion
max_iter = 10 # Maximal number of iterations between updating time shifts and intensities in the centering step
t_vec = results[[1]]$network_param$t_vec

### Apply proposed initialization scheme
res = get_init(edge_time_mat_list = edge_time_mat_list,
               N_clus = N_clus_tmp,
               t_vec = t_vec)

clusters_list_init = res$clusters_list # Initial clustering
n0_vec_list_init = res$n0_vec_list # Initial node-specific time shifts
n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init) # Initial edge-specific time shifts

res = do_cluster_cdf(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                     total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                     clusters_list_init = clusters_list_init,
                     n0_vec_list_init = n0_vec_list_init, 
                     n0_mat_list_init = n0_mat_list_init,
                     freq_trun = freq_trun,
                     conv_thres = conv_thres,
                     MaxIter = MaxIter)

set.seed(183)
mem = clus2mem(clusters_list_init[[1]])
mem2 = mem
mem[sample(1:30,6,replace = FALSE)] = sample(c(1:3), 6, replace=TRUE)
clusters_list_init_2 = list(mem2clus(mem))

n0_vec_list_init_2 = n0_vec_list_init
n0_vec_list_init_2[[1]] = round(jitter(n0_vec_list_init[[1]]))
n0_mat_list_init_2 = n0_vec2mat(n0_vec = n0_vec_list_init_2) 

### Apply algorithm 
res = do_cluster_cdf(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                     total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                     clusters_list_init = clusters_list_init_2,
                     n0_vec_list_init = n0_vec_list_init_2, 
                     n0_mat_list_init = n0_mat_list_init_2,
                     freq_trun = freq_trun,
                     conv_thres = conv_thres,
                     MaxIter = MaxIter)

clusters_list_est = res$clusters_list # Estimated clusters
v_vec_list_est = res$v_vec_list # Estimated time shifts
center_cdf_array_est = res$center_cdf_array # Estimated cumulative connecting intensities
center_pdf_array_est = res$center_pdf_array # Estimated connecting intensities


cdf_true_array = network_list$cdf_true_array # True cumulative connecting intensities


### Match the estimated clusters and the true clusters 
### Identifiable upon to permutation
### ONLY for simulation 

permn = find_permn(center_cdf_array_from = center_cdf_array_est, 
                   center_cdf_array_to = cdf_true_array)$permn


### Plot estimated intensities (black curves) and true intensities (red curves)

center_cdf_array_est_permn = center_cdf_array_est[permn, permn, ] # Estimated cumulative intensities with permuted cluster label

plot_pdf_array(pdf_array_list = center_cdf_array_est_permn, 
               pdf_true_array = cdf_true_array, 
               t_vec = t_vec, 
               y_lim = c(0,max(c(center_cdf_array_est_permn, cdf_true_array))) )

