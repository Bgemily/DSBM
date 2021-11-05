# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Load simulation result and get network parameters -----------------------

load("../Results/Rdata/SNR_Vis0/apply_ppsbm_ICL/pr=0.4,n=30,beta=1.3/n/90/N_trial5_20211028_203955.Rdata")
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
N_clus_min=N_clus-2 
N_clus_max=N_clus+2

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


# Apply PPSBM -------------------------------------------------------------

library(ppsbm)
edge_time_mat = edge_time_mat_list[[1]]
time.seq = numeric(sum(edge_time_mat<Inf))
type.seq = numeric(sum(edge_time_mat<Inf))
current_ind = 1
for (i in 1:nrow(edge_time_mat)) {
  for (j in i:ncol(edge_time_mat)) {
    if (edge_time_mat[i,j]<Inf){
      time.seq[current_ind] = edge_time_mat[i,j]
      type.seq[current_ind] = convertNodePair(i, j, n = nrow(edge_time_mat), directed = FALSE)
      current_ind = current_ind+1
    }
  }
}
data = list(time.seq=time.seq, type.seq=type.seq, Time=total_time)
Nijk = statistics(data, nrow(edge_time_mat), K=2^6, directed = FALSE)

res_list_ppsbm = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), d_part=5, 
                         Qmin=N_clus_min, Qmax=N_clus_max, directed=FALSE, sparse=FALSE, method="hist")

save(res_list_ppsbm, network_param,
     file = '../Results/Rdata/SNR_Vis0/apply_ppsbm_ICL/pr=0.4,n=90,beta=1.3,one_instance.Rdata')


# Get ppsbm ICL -----------------------------------------------------------

# ICL-model selection
data = list(Nijk=Nijk, Time=total_time)
sol.selec_Q <- modelSelection_Q(data=data,
                                n=nrow(edge_time_mat),
                                Qmin=N_clus_min,
                                Qmax=N_clus_max,
                                directed=FALSE,
                                sparse=FALSE,
                                sol.hist.sauv=res_list_ppsbm)

# best number Q of clusters:
loglik_vec_ppsbm = sol.selec_Q$all.compl.log.likelihood

# Apply our method --------------------------------------------------------

res_list = list()
for (N_clus_tmp in N_clus_min:N_clus_max) {
  ### Get initialization -----------
  res = get_init_v3(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp, 
                    v_true_list = v_true_list, 
                    jitter_time_rad = 0,
                    t_vec = t_vec)
  
  clusters_list_init = res$clusters_list
  n0_vec_list_init = res$n0_vec_list
  n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
  
  # Apply algorithm ---------
  
  
  ### Estimation z,v,f based on cdf
  res = do_cluster_v8.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                        total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                        clusters_list_init = clusters_list_init,
                        n0_vec_list_init = n0_vec_list_init, n0_mat_list_init = n0_mat_list_init,
                        )
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$n0_vec_list -> n0_vec_list_est
  res$n0_mat_list -> n0_mat_list_est
  res$center_cdf_array -> center_cdf_array_est
  
  
  
  ### Estimation z,v,f based on pdf
  res = do_cluster_v14.2.1(edge_time_mat_list = edge_time_mat_list, N_clus = N_clus_tmp,
                           clusters_list_init = clusters_list_est,
                           n0_vec_list_init = n0_vec_list_est, n0_mat_list_init = n0_mat_list_est,
                           total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                           freq_trun=freq_trun, 
                           conv_thres=conv_thres, MaxIter=MaxIter,
                           opt_radius=opt_radius,
                           )
  
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$loss_history -> loss_history
  res$align_time -> align_time
  res$cluster_time -> cluster_time
  res$center_pdf_array -> center_pdf_array_est
  
  # Save results of N_clus_tmp ----------------------------------------------
  
  res_list = c(res_list, list(res))
  
}

save(res_list, network_param,
     file = '../Results/Rdata/SNR_Vis0/main_v5_v7_largefreqtrun/pr=0.4,n=90,beta=1.3,one_instance.Rdata')

# Select best cluster number using ICL ------------------------------------

sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list, 
                           N_node_vec = network_param$N_node_vec, 
                           N_clus_min = N_clus_min, 
                           N_clus_max = N_clus_max, 
                           result_list = res_list, 
                           total_time = total_time)

N_clus_est = sel_mod_res$N_clus_est
ICL_vec = sel_mod_res$ICL_vec 
compl_log_lik_vec = sel_mod_res$compl_log_lik_vec 
penalty_vec = sel_mod_res$penalty_vec


# Check clustering result -------------------------------------------------

for (tmp in 1:3) {
  print(paste0("When the number of clusters = ", tmp," , the estimated clusters are: "))
  clusters = mem2clus(apply(res_list_ppsbm[[tmp]]$tau, 2, which.max), N_clus_min = tmp) 
  print(clusters)
}

for (tmp in 1:3) {
  print(paste0("When the number of clusters = ", tmp," , the estimated clusters are: "))
  clusters = res_list[[tmp]]$clusters_list[[1]]
  print(clusters)
}


# Plot estimated intensities ----------------------------------------------

plot_pdf_array_v2(pdf_array_list = list(res_list[[3]]$center_pdf_array), 
                  y_lim = c(0,0.04),
                  t_vec = t_vec)

plot(x=seq(0,200,length.out=ncol(res_list_ppsbm[[3]]$logintensities.ql)), 
     y=exp(res_list_ppsbm[[1]]$logintensities.ql[1,])*I(res_list_ppsbm[[1]]$logintensities.ql[1,]!=0), 
     type='s',col=2)

lines(x=t_vec, y=res_list[[1]]$center_pdf_array[1,1,],type='l')

