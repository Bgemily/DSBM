

### Calculate the gradient of loss function 
### [WARNING: the gradient might be WRONG when there are multiple subjects!]
get_gradient = function(edge_time_mat_list, clusters_list, 
                        n0_vec_list, n0_mat_list=NULL,
                        freq_trun=5, t_vec=seq(0, 200, length.out=1000) )
{
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_node_vec = sapply(edge_time_mat_list, nrow)
  N_clus = length(clusters_list[[1]])
  membership_list = lapply(clusters_list, clus2mem)
  
  
  if (is.null(n0_mat_list)) {
    n0_mat_list = n0_vec2mat(n0_vec = n0_vec_list)
  }
  
  
  ### Compute center_fft_array, i.e. phi_qk's
  center_fft_array = get_center_fft_array(edge_time_mat_list = edge_time_mat_list, 
                                          clusters_list = clusters_list, 
                                          freq_trun = freq_trun, 
                                          n0_mat_list = n0_mat_list, 
                                          rmv_conn_prob = TRUE,
                                          t_vec = t_vec)
  
  
  ### Compute gradient for all time shifts v_{m,i}'s
  gradient_vec_list = lapply(N_node_vec, numeric)
  for (m in 1:N_subj) {
    for (i in N_node_vec[m]:1) {
      
      q = membership_list[[m]][i] #z_mi
      A = matrix(nrow=N_clus, ncol=2*freq_trun+1)
      B = matrix(nrow=N_clus, ncol=2*freq_trun+1)
      for (k in 1:N_clus) {
        ### Define frequency vector ###
        l_vec = 0:(2*freq_trun)
        l_vec = c(tail(l_vec, freq_trun)-length(l_vec),
                  head(l_vec, freq_trun+1))
        ### ### ###
        
        ### Compute A[k,] ###
        J_mki = which(membership_list[[m]]==k & 
                        n0_vec_list[[m]]<=n0_vec_list[[m]][i] & 
                        edge_time_mat_list[[m]][i,]<Inf)
        J_mki = setdiff(J_mki, i)
        
        tmp = edge_time_mat_list[[m]][i,J_mki]
        if (length(tmp)==0) {
          c_bar_mki_vec = rep(0,2*freq_trun+1)
        }
        else{
          ### Get empirical intensity of event times
          event_time_vec = tmp
          emp_intens_vec = hist(event_time_vec, breaks=t_vec, plot=FALSE)$counts / length(event_time_vec)
          emp_intens_vec = c(0,emp_intens_vec)
          ### Get normalized fourier series
          emp_intens_fft = fft(emp_intens_vec) / length(t_vec)
          ### Get truncated fourier series
          emp_intens_fft_trun = c(tail(emp_intens_fft, freq_trun),
                                  head(emp_intens_fft, freq_trun+1))
          c_bar_mki_vec = emp_intens_fft_trun
        }
        
        
        A[k,] = 1i*2*pi*(l_vec/length(t_vec)) * length(J_mki) *
          Conj( c_bar_mki_vec * exp(1i*2*pi*n0_vec_list[[m]][i]*(l_vec/length(t_vec))) ) 
        ### ###
        
        ### Compute B[k,] ###
        B[k,] = center_fft_array[q,k,] 
        ### ###
        
        
        
      } 
      
      gradient = length(t_vec) * sum( 2 * Re(A*B) ) ### "length(t_vec)" comes from \|f\|^2 = T * \sum |c_l|^2
      gradient_vec_list[[m]][i] = gradient
      
      
    }
  }
  
  
  return(gradient_vec_list)
  
  
}

