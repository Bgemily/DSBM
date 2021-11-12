
# -------------------------------------------------------------------------

# Move f1 towards f2 by n0. Positive n0: towards right. Negative n0: towards left.

align_multi_fft_gd = function(fft_origin_list, fft_target_list, 
                              n0=0, step_size=0.02, shrink = 0.5,
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
  theta_prime_list = fft_origin_list
  gamma_prime_list = fft_target_list
  
  
  ### Gradient descent
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  n0_curr = n0_upd = n0
  while (dist_redu>stopping_redu && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd_list = mapply(FUN = gradient_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, 
                     MoreArgs = list(n0=n0_curr, pad=0), SIMPLIFY = FALSE)
    if(is.null(weights)) 
      gd = mean(unlist(gd_list))
    else
      gd = sum(unlist(gd_list)*weights)
    
    ### Set initial learning rate
    step_size = (length(t_vec)/4) / gd

    ### Update n0
    n0_upd = n0_curr - (step_size)*gd
    n0_upd = round(n0_upd) 
    if (n0_upd<n0_min) {
      n0_upd = n0_min
    }
    if (n0_upd>n0_max) {
      n0_upd = n0_max
    }
    
    ### Evaluate loss (dist_upd)
    dist_list = mapply(FUN = distance_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, 
                       MoreArgs = list(n0 = n0_upd, pad=0), SIMPLIFY = FALSE)
    if(is.null(weights))
      dist_upd = sqrt(mean(unlist(dist_list)^2))
    else
      dist_upd = sqrt(sum(weights*(unlist(dist_list))^2)) 
    
    ### Adjust step size
    # TODO: finish this step
    N_shrink = 0
    while (dist_upd>dist_curr & N_shrink<=10) {
      ### Shrink step_size
      N_shrink = N_shrink+1
      step_size = step_size*shrink
      
      ### Estimate n0 based on shrunk step size
      n0_upd = n0_curr - (step_size)*gd
      n0_upd = round(n0_upd) 
      if (n0_upd<n0_min) {
        n0_upd = n0_min
      }
      if (n0_upd>n0_max) {
        n0_upd = n0_max
      }
      
      ### Evaluate loss (dist_upd)
      dist_list = mapply(FUN = distance_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, 
                         MoreArgs = list(n0 = n0_upd, pad=0), SIMPLIFY = FALSE)
      if(is.null(weights))
        dist_upd = sqrt(mean(unlist(dist_list)^2))
      else
        dist_upd = sqrt(sum(weights*(unlist(dist_list))^2)) 
      
    }
    
    if (N_shrink>10) {
      message("[align_multi_fft_gd]: Reached maximum shrinkage number.")
    }
    
    ### Compute reduction in loss(distance)
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    
    # browser()
    # print("n0_curr:");print(n0_curr)
    # print("n0_upd:");print(n0_upd)
    # print("(step_size)*gd");print((step_size)*gd)
    
    ### *upd -> *curr
    dist_curr = dist_upd
    n0_curr = n0_upd
    
  }
  
  
  if (iter_count==MaxIter) {
    warning("[align_multi_fft_gd_v2]: Reached max iteration number. Consider adjusting the step size.")
  }
  
  return(list(n0=n0_upd, dist = dist_curr))
  
  
}





# Test --------------------------------------------------------------------


# t_vec=seq(0,100,0.1)
# freq_trun=5
# 
# samp1 = samp2 = rep(0,length(t_vec))
# samp1[100:300+0] = 1
# samp2[100:300+400] = 1
# 
# plot(samp1,type='l')
# plot(samp2,type='l')
# 
# fft1 = fft(samp1)/length(t_vec)
# fft2 = fft(samp2)/length(t_vec)
# 
# fft1_trun = c(head(fft1, freq_trun+1),
#          rep(0,length(fft1)-2*freq_trun-1),
#          tail(fft1, freq_trun))
# fft2_trun = c(head(fft2, freq_trun+1),
#          rep(0,length(fft2)-2*freq_trun-1),
#          tail(fft2, freq_trun))
# f1 = Re(fft( fft1_trun , inverse = TRUE))
# f2 = Re(fft( fft2_trun , inverse = TRUE))
# plot(f1, type='l')
# lines(f2, col=2)
# 
# 
# res = align_multi_fft_gd(fft_origin_list = list(fft1_trun), fft_target_list = list(fft2_trun),
#                          n0 = 50, t_vec = t_vec, freq_trun = freq_trun, step_size = 100)
# res$n0
# 
# res = align_multi_curves_gd_v2(f1,f2,n0_min = 0,n0 = 200)
# res$n0
# 
