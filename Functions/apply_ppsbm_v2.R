
### Use generate_network
### Output measurements of errors
apply_ppsbm_v2 = function(### Parameters for generative model
                          SEED, N_subj=1, N_node_vec = rep(90,N_subj),
                          N_clus=3, clus_size_mat = matrix(N_node_vec/N_clus, nrow=N_subj, ncol=N_clus),
                          total_time=200, 
                          conn_patt_var=1, conn_patt_sep = 1.5, const=20, conn_prob_mean = 1, conn_prob_rad = 0, 
                          time_shift_struc=max, time_shift_mean_vec = rep(20,N_clus), 
                          time_shift_rad = min(time_shift_mean_vec),
                          ### Parameters for algorithms
                          Qmin=N_clus, Qmax=N_clus, bw=5,
                          t_vec = seq(0,total_time,length.out=1000),
                          ### Result saving setting
                          save_center_pdf_array=FALSE
                          )
{
  
  # Generate networks -------------------------------------------------------
  
  
  network_list = generate_network(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                                      N_clus = N_clus, clus_size_mat = clus_size_mat,
                                      total_time = total_time, 
                                      conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                                      conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                                      time_shift_struc = time_shift_struc,
                                      t_vec = t_vec,
                                      time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  
  edge_time_mat_list = network_list$edge_time_mat_list
  cdf_true_array = network_list$cdf_true_array
  pdf_true_array = network_list$pdf_true_array
  membership_true_list = network_list$membership_true_list
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$time_shift_list
  order_true_list = lapply(v_true_list, function(v_vec)order(v_vec))
  
  
  edge_time_mat = edge_time_mat_list[[1]]
  
  

  # Apply PPSBM -------------------------------------------------------------

  library(ppsbm)
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
  
  time_start = Sys.time()
  # res = mainVEM(data=data, n=nrow(edge_time_mat), Qmin=3, directed=FALSE, method="kernel",
  #               d_part=0, n_perturb=0)
  res = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), d_part=5, 
                Qmin=Qmin, Qmax=Qmax, directed=FALSE, method="hist")[[1]]
  time_end = Sys.time()
  time_estimation = time_end - time_start
  time_estimation = as.numeric(time_estimation, units='secs')
  res$clusters = mem2clus(apply(res$tau, 2, which.max)) 
  
  
  clusters_list_est = list(res$clusters)
  
  ### Extract estimated intensities 
    N_clus_est = length(clusters_list_est[[1]])
    center_pdf_array_est = array(dim = c(N_clus_est,N_clus_est,length(t_vec)))
    ind_qk = 1
    for (q in 1:N_clus_est) {
      for (k in q:N_clus_est) {
        ### Extract estimated intensity for cluster pair (q,k)
        intensity_tmp = exp(res$logintensities.ql[ind_qk, ])
        ### Fix bug for ppsbm
        intensity_tmp[intensity_tmp==1] = 0
        ### Add breakpoints in time grid
        rep_time = floor( (1:length(intensity_tmp))*length(t_vec)/length(intensity_tmp) ) -
                    floor( (0:(length(intensity_tmp)-1))*length(t_vec)/length(intensity_tmp) )
        intensity_tmp = rep(intensity_tmp, time = rep_time)

        center_pdf_array_est[q,k, ] = intensity_tmp
        center_pdf_array_est[k,q, ] = center_pdf_array_est[q,k, ]
        ind_qk = ind_qk + 1
      }
    }
  
  ### Estimate intensities using kernel smoothing
  center_pdf_array_est_kernel = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list, 
                                                 clusters_list = clusters_list_est, 
                                                 n0_mat_list = list(matrix(0,nrow(edge_time_mat), 
                                                                           ncol(edge_time_mat))), 
                                                 t_vec = t_vec)
  
  
  # Compute errors of clusters, i.e. Z ------------------------------------
  
  ARI_vec = numeric(length = N_subj)
  for (m in 1:N_subj) {
    ARI_tmp = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est[[m]]), 
                          memb_true_vec = membership_true_list[[m]])
    ARI_vec[m] = ARI_tmp
  }
  ARI_mean = mean(ARI_vec)
  
  
  # Compute errors of conn patts, i.e. F ------------------------------------
  
  ### Match clusters
  membership_est_vec = clus2mem(clusters_list_est[[1]])
  z_hat = t(dummies::dummy(membership_est_vec))
  z_true = t(dummies::dummy(membership_true_list[[1]]))
  permn = ppsbm::permuteZEst(z = z_true, hat.z = z_hat)
  if(nrow(z_hat)<nrow(z_true)){
    permn_tmp = rep(setdiff(1:nrow(z_true),permn),nrow(z_true))
    permn_tmp[sort(unique(apply(res$tau, 2, which.max)))] = permn
    permn = permn_tmp
  }
  
  ### Calculate imse for center_pdf_array_est 
  center_pdf_array_est_permn = center_pdf_array_est[permn, permn, ]
  dist_mat = matrix(nrow=N_clus, ncol=N_clus)
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      dist_mat[q,k] = sqrt(sum( (center_pdf_array_est_permn[q,k,] - pdf_true_array[q,k,])^2 * 
                                  (t_vec[2]-t_vec[1]) ))
    }
  }
  F_mean_sq_err = mean(dist_mat[upper.tri(dist_mat, diag=TRUE)]^2)
  
  ### Calculate imse for center_pdf_array_est_kernel
  center_pdf_array_est_kernel_permn = center_pdf_array_est_kernel[permn, permn, ]
  dist_mat_kernel = matrix(nrow=N_clus, ncol=N_clus)
  for (q in 1:N_clus) {
    for (k in 1:N_clus) {
      dist_mat_kernel[q,k] = sqrt(sum( (center_pdf_array_est_kernel_permn[q,k,] - pdf_true_array[q,k,])^2 * 
                                  (t_vec[2]-t_vec[1]) ))
    }
  }
  F_mean_sq_err_kernel = mean(dist_mat_kernel[upper.tri(dist_mat_kernel, diag=TRUE)]^2)
  
  
  # Extract network related parameters -----------------------------------------
  
  network_param = list(SEED = SEED, N_subj = N_subj, N_node_vec = N_node_vec, 
                       N_clus = N_clus, clus_size_mat = clus_size_mat,
                       total_time = total_time, t_vec = t_vec,
                       conn_patt_var = conn_patt_var, conn_patt_sep = conn_patt_sep, const = const,
                       conn_prob_mean = conn_prob_mean, conn_prob_rad = conn_prob_rad, 
                       time_shift_struc = time_shift_struc,
                       time_shift_mean_vec = time_shift_mean_vec, time_shift_rad = time_shift_rad)
  
  
  # Output ------------------------------------------------------------------
  
  return(list(network_param=network_param, 
              t_vec=t_vec,
              time_estimation=time_estimation,
              clusters_list_est=clusters_list_est,
              ARI_vec=ARI_vec, ARI_mean=ARI_mean,
              F_mean_sq_err_kernel=F_mean_sq_err_kernel,
              F_mean_sq_err=F_mean_sq_err,
              center_pdf_array_est=switch(save_center_pdf_array, 
                                          yes=center_pdf_array_est,
                                          no=NULL)
              ))
}


