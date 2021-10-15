rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(ppsbm)

ppsbm_vec = c()
our_vec = c()
SEED_vec = c()

for (trial in 1:50) {
  
  SEED=sample(1e4,1)
  SEED_vec = c(SEED_vec,SEED)
  total_time=200
  t_vec = seq(0,200,length.out=200)
  network_list = generate_network2_v3(SEED = SEED, N_subj = 1, N_node_vec = c(30), 
                                      N_clus = 3, conn_prob_mean = 0.4, 
                                      conn_patt_sep = 1.8,
                                      t_vec = seq(0,200,length.out=200),
                                      total_time = 200,
                                      time_shift_mean_vec = rep(0,3))
  
  edge_time_mat_list = network_list$edge_time_mat_list
  edge_time_mat = edge_time_mat_list[[1]]
  
  # Apply PPSBM -------------------------------------------------------------
  
  time.seq = numeric(sum(edge_time_mat<Inf))
  type.seq = numeric(sum(edge_time_mat<Inf))
  current_ind = 1
  for (i in 1:nrow(edge_time_mat)) {
    for (j in 1:ncol(edge_time_mat)) {
      if (edge_time_mat[i,j]<Inf){
        time.seq[current_ind] = edge_time_mat[i,j]
        type.seq[current_ind] = convertNodePair(i, j, n = nrow(edge_time_mat), directed = FALSE)
        current_ind = current_ind+1
      }
    }
  }
  data = list(time.seq=time.seq, type.seq=type.seq, Time=total_time)
  Nijk = statistics(data, nrow(edge_time_mat), K=2^6, directed = FALSE)
  
  res = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), d_part=5, 
                Qmin=1, Qmax=1, directed=FALSE, sparse=FALSE, method="hist")
  
  # ICL-model selection
  data = list(Nijk=Nijk, Time=total_time)
  sol.selec_Q <- modelSelection_Q(data=data,
                                  n=nrow(edge_time_mat),
                                  Qmin=1,
                                  Qmax=1,
                                  directed=FALSE,
                                  sparse=FALSE,
                                  sol.hist.sauv=res)
  log_lik_vec = sol.selec_Q$all.compl.log.likelihood
  
  
  
  
  # Calculate log_lik_vec using our method --------------------------------------------------------
  
  res_reformat = list()
  
  res_reformat$clusters_list = list( mem2clus(apply(res[[1]]$tau, 2, which.max),N_clus_min = 1) )
  res_reformat$v_vec_list = list( rep(0,nrow(edge_time_mat)) )
  
  ### Extract estimated intensities
  N_clus_est = 1
  center_pdf_array_est = array(dim = c(N_clus_est,N_clus_est,length(t_vec)))
  ind_qk = 1
  for (q in 1:N_clus_est) {
    for (k in q:N_clus_est) {
      ### Extract estimated intensity for cluster pair (q,k)
      intensity_tmp = exp(res[[1]]$logintensities.ql[ind_qk, ])
      ### Add breakpoints in time grid
      rep_time = floor( (1:length(intensity_tmp))*length(t_vec)/length(intensity_tmp) ) - 
        floor( (0:(length(intensity_tmp)-1))*length(t_vec)/length(intensity_tmp) )
      intensity_tmp = rep(intensity_tmp, time = rep_time)
      
      center_pdf_array_est[q,k, ] = intensity_tmp
      center_pdf_array_est[k,q, ] = center_pdf_array_est[q,k, ]
      ind_qk = ind_qk + 1
    }
  }
  res_reformat$center_pdf_array = center_pdf_array_est
  
  
  
  res_list = list(res_reformat)
  sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list, 
                             N_node_vec = c(30), 
                             N_clus_min = 1, 
                             N_clus_max = 1, 
                             result_list = res_list, 
                             t_vec = seq(0,200,length.out=200), 
                             freq_trun = 10)
  
  compl_log_lik_vec = sel_mod_res$compl_log_lik_vec 
  
  
  ppsbm_vec = c(ppsbm_vec, log_lik_vec)
  our_vec = c(our_vec, compl_log_lik_vec)
  
  
}

pdf(file="../Results/Plots/Temp/log_lik_our_vs_ppsbm.pdf", 
    width = 4, height = 4)
plot(ppsbm_vec,our_vec)
dev.off()

mean(our_vec) - mean(ppsbm_vec)


