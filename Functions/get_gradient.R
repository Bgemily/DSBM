

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
        J_mki = which(membership_list[[m]]==k & n0_vec_list[[m]]<=n0_vec_list[[m]][i])
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



# Test --------------------------------------------------------------------

# edge_time_mat = kronecker(matrix(1,2,2), matrix(100,5,5))
# t_vec = seq(0,200,1)
# t_unit = t_vec[2] - t_vec[1]
# n0_mat = (edge_time_mat*0+5)/t_unit
# n0_vec = rep(5,nrow(edge_time_mat))
# edge_time_mat = edge_time_mat + 5
# 
# freq_trun = 5
# tmp = 0*t_vec
# tmp[100] = 1
# fft_true_array = fft(tmp)/length(t_vec)
# fft_true_array = c(tail(fft_true_array, freq_trun),
#                    head(fft_true_array, freq_trun+1))
# 
# clusters = list(1:5,6:10)
# 
# edge_time_mat_list = list(edge_time_mat,edge_time_mat)
# n0_vec_list = list(n0_vec, n0_vec)
# n0_mat_list = list(n0_mat, n0_mat)
# clusters_list = list(clusters,clusters)
# center_fft_array = array(dim=c(2,2,2*freq_trun+1))
# center_fft_array[1,1,] = center_fft_array[2,2,] = center_fft_array[1,2,] = center_fft_array[2,1,] = fft_true_array
# 
# 
# get_gradient(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list, 
#              n0_vec_list = lapply(n0_vec_list, function(x)x+runif(length(x),0,10)), 
#              freq_trun = freq_trun, t_vec = t_vec)

