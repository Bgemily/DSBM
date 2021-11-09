# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Load simulation result and get network parameters -----------------------

load("../Results/Rdata/SNR_Vnot0/main_v5_pdf_v2/freqtrun/7/pr=1,n=30,beta=1.2/n/90/N_trial5_20211108_154619.Rdata")
network_param = results[[1]]$network_param

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

freq_trun=7 
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
  
  clusters_list_init -> clusters_list_est
  n0_vec_list_init -> n0_vec_list_est
  n0_mat_list_init -> n0_mat_list_est
  
  ### Estimation z,v,f based on pdf
  res = do_cluster_v14.2.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                           clusters_list_init = clusters_list_est,
                           n0_vec_list_init = n0_vec_list_est, n0_mat_list_init = n0_mat_list_est,
                           total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                           freq_trun=freq_trun, 
                           conv_thres=conv_thres, MaxIter=MaxIter,
                           opt_radius=opt_radius)
  
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$loss_history -> loss_history
  res$align_time -> align_time
  res$cluster_time -> cluster_time
  res$center_pdf_array -> center_pdf_array_est
  
  
  ### Get estimated pdf using kernel smoothing
  v_mat_list_est = n0_vec2mat(n0_vec = v_vec_list_est)
  n0_mat_list_est = lapply(v_mat_list_est, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
  center_pdf_array_est = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list,
                                                 clusters_list = clusters_list_est,
                                                 n0_mat_list = n0_mat_list_est,
                                                 t_vec = t_vec)
  res$center_pdf_array = center_pdf_array_est
  
  
  # Save results of N_clus_tmp ----------------------------------------------
  
  res_list = c(res_list, list(res))
  
}

# save(res_list, network_param,
#      file = '../Results/Rdata/SNR_Vis0/main_v5_v7_largefreqtrun/pr=0.4,n=90,beta=1.3,one_instance.Rdata')



# Plot estimated time shifts ----------------------------------------------

plot(y=res$v_vec_list[[1]],x=v_true_list[[1]])
abline(a=0,b=1,col='red')

# Plot estimated intensities ----------------------------------------------

plot_pdf_array_v2(pdf_array_list = list(res_list[[1]]$center_pdf_array), 
                  pdf_true_array = network_list$pdf_true_array,
                  clus_size_vec = sapply(res$clusters_list[[1]],length),
                  y_lim = c(0,0.06),
                  t_vec = t_vec)

plot(x=seq(0,200,length.out=ncol(res_list_ppsbm[[3]]$logintensities.ql)), 
     y=exp(res_list_ppsbm[[1]]$logintensities.ql[1,])*I(res_list_ppsbm[[1]]$logintensities.ql[1,]!=0), 
     type='s',col=2)

lines(x=t_vec, y=res_list[[1]]$center_pdf_array[1,1,],type='l')

