# Move f1 towards f2 by n0. Positive n0: towards right. Negative n0: towards left.


align_multi_fft_gd = function(fft_origin_list, fft_target_list, 
                              n0=0, step_size=0.02, 
                              freq_trun = 10,
                              MaxIter=1000, stopping_redu=0.0001, 
                              weights=NULL, t_vec=seq(0,200,length.out=1000),
                              n0_min = 0, n0_max = length(t_vec))
{
  if(!is.list(fft_origin_list)) fft_origin_list = list(fft_origin_list)
  if(!is.list(fft_target_list)) fft_target_list = list(fft_target_list)
  
  ### Prevent n0_min and n0_max from further changing
  n0_min -> n0_min
  n0_max -> n0_max 
  
  
  ### Check weights' validity
  if(!is.null(weights)&&sum(weights)==0) stop("Weights cannot sum up to zero.")
  if(!is.null(weights)&&sum(weights)!=1) weights = weights/sum(weights)
  
  
  ### Compute terms needed in gradients
  theta_prime_list = lapply(fft_origin_list, function(fft)c(tail(fft,freq_trun+1), 
                                                            head(fft,freq_trun))) 
  gamma_prime_list = lapply(fft_target_list, function(fft)c(tail(fft,freq_trun+1), 
                                                            head(fft,freq_trun)))
  
  
  ### Gradient descent
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  while (dist_redu>stopping_redu && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd_list = mapply(FUN = gradient_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, 
                     MoreArgs = list(n0=n0, pad=0), SIMPLIFY = FALSE)
    if(is.null(weights)) 
      gd = mean(unlist(gd_list))
    else
      gd = sum(unlist(gd_list)*weights)
    
    n0 = n0 - (step_size)*gd
    n0 = round(n0) 
    
    print((step_size)*gd)
    
    if (n0<n0_min) {
      n0 = n0_min
    }
    if (n0>n0_max) {
      n0 = n0_max
    }
    
    dist_list = mapply(FUN = distance_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, 
                       MoreArgs = list(n0 = n0, pad=0), SIMPLIFY = FALSE)
    if(is.null(weights))
      dist_upd = sqrt(mean(unlist(dist_list)^2))
    else
      dist_upd = sqrt(sum(weights*(unlist(dist_list))^2)) 
    
    
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    dist_curr = dist_upd
  }
  
  
  if (iter_count==MaxIter) {
    warning("[align_multi_fft_gd_v2]: Reached max iteration number. Consider adjusting the step size.")
  }
  
  return(list(n0=n0, dist = dist_curr))
  
  
}





# Test --------------------------------------------------------------------


# t_vec=seq(0,200,length.out=1000)
# freq_trun=5
# 
# samp1 = samp2 = rep(0,length(t_vec))
# samp1[100:300+0] = 1
# samp2[100:300+200] = 1
# 
# 
# fft1 = fft(samp1)/length(t_vec)
# fft2 = fft(samp2)/length(t_vec)
# 
# fft1 = c(tail(fft1, freq_trun),
#          head(fft1, freq_trun+1))
# fft2 = c(tail(fft2, freq_trun),
#          head(fft2, freq_trun+1))
# 
# plot(Re(fft( c(tail(fft1, freq_trun+1),
#                rep(0, length(t_vec)-2*freq_trun-1),
#                head(fft1, freq_trun)) , inverse = TRUE)))
# 
# lines(Re(fft( c(tail(fft2, freq_trun+1),
#                rep(0, length(t_vec)-2*freq_trun-1),
#                head(fft2, freq_trun)) , inverse = TRUE)), col=2)
# 
# 
# res = align_multi_fft_gd(fft_origin_list = list(fft1), fft_target_list = list(fft2),
#                          n0 = 200, t_vec = t_vec, freq_trun = freq_trun, step_size = 100)
# res$n0
# res$dist

