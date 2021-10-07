
# Shift curve -----------------------------------------------

### n0>0: towards right. n0<0: towards left.
shift_v2 = function(f_origin, n0, pad = tail(f_origin,1), pp=FALSE)
{
  n0 = round(n0)
  N = length(f_origin)
  if (!pp & n0<=0)
    return(c(f_origin[(abs(n0)+1):N], rep(pad, abs(n0))))
  else if (n0<=0)
    return(c(f_origin[(abs(n0)+1):N], rep(0, abs(n0))))
  else
    return( c(rep(0, abs(n0)), head(f_origin, N-abs(n0))) )
}

# ### Test
# tmp = c(1:100*0,1:100,1:100*0+100)
# plot(tmp)
# plot(shift_v2(tmp,-50))


# Compute distance with given shift ---------------------------------------

distance_v2 = function(theta_prime, gamma_prime, n0, t_unit=0.05, pad=1) # L2 distance
{
  N = length(theta_prime)
  l_vec = seq(1,N-1)
  
  d = sum(abs(theta_prime[l_vec+1]*exp(1i*2*pi*l_vec*(-n0)/N)-gamma_prime[l_vec+1])^2) + 
        abs(theta_prime[1]+(-n0)*pad-gamma_prime[1])^2
  d = sqrt(t_unit / N * d)
  return(d)
}



# gradient -----------------------------------------------------------------

# loss function = squared L2 distance*(N/t_unit)
gradient_v2 = function(theta_prime, gamma_prime, n0, pad=1)
{
  N = length(theta_prime)
  if(N != length(gamma_prime)) stop("Length of theta and gamma do not match.")
  
  theta_prime_0 = theta_prime[1]
  gamma_prime_0 = gamma_prime[1]
  
  theta_prime = theta_prime[-1]
  gamma_prime = gamma_prime[-1]
  theta_prime = c(tail(theta_prime, length(theta_prime)/2), 
                  head(theta_prime, length(theta_prime)/2))
  gamma_prime = c(tail(gamma_prime, length(gamma_prime)/2), 
                  head(gamma_prime, length(gamma_prime)/2))
  
  l_vec = 1:(N-1)
  l_vec = c(tail(l_vec, length(l_vec)/2)-N,
            head(l_vec, length(l_vec)/2))
  
  gradient = 2 * sum(Re(1i*2*pi*(l_vec/N)*(theta_prime)*Conj(gamma_prime)*exp(1i*2*pi*(-n0)*l_vec/N))) + 
              (-2)*Re(theta_prime_0 - n0*pad - gamma_prime_0)*pad 
    
  return (gradient)
}

# grads = sapply(seq(0,400,0.1),function(x)gradient(theta_prime, gamma_prime, x))
# plot(seq(0,400,0.1), grads, type = 'l')
# abline(h=0,col=2)


# Align curve (using gradient) -------------------------------------------------------------

# Non-negative n0. Move f_origin towards f_shift by n0. If pp=FALSE, pad f_origin[length(f_origin)] at the end of f_origin when shifting; otherwise, pad 0.
align_curves_gd_ = function(f_origin, f_shift, n0, step_size, MaxIter=1000, stopping_redu=0.01, pp=FALSE, t_unit=0.05)
{
  ########### V1
  # theta = fft(f_origin)
  # gamma = fft(f_shift)
  # N = length(theta)
  # k = seq(1, N-1)
  # pad = f_origin[length(f_origin)]
  # {
  #   if (!pp){
  #     ############# V1
  #     # theta_prime = c(theta[1], theta[2:N]+1/(1-exp(-1i*2*pi*k/N)))
  #     # gamma_prime = c(gamma[1], gamma[2:N]+1/(1-exp(-1i*2*pi*k/N)))
  #     ############# V2
  #     theta_prime = c(theta[1], theta[2:N]+1*pad/(1-exp(-1i*2*pi*k/N)))
  #     gamma_prime = c(gamma[1], gamma[2:N]+1*pad/(1-exp(-1i*2*pi*k/N)))
  #     ########################################
  #   }
  #   else{
  #     theta_prime = theta
  #     gamma_prime = gamma
  #   }
  # }
  ############ V2
  res = get_theta_gamma_prime(f_origin = f_origin, f_shift = f_shift, pp = pp)
  theta_prime = res$theta_prime
  gamma_prime = res$gamma_prime
  ########################
  
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  while (dist_redu>stopping_redu && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd = gradient(theta_prime = theta_prime, gamma_prime = gamma_prime, n0 = n0, pp=pp, pad = pad); 
    n0 = n0 - step_size*gd
    n0 = round(n0) 
    
    if(n0<0) { 
      n0 = 0; 
      dist_curr = distance(theta_prime = theta_prime, gamma_prime = gamma_prime, n0 = n0, pp=pp, t_unit = t_unit, pad = pad)
      break 
    }
    if(n0>length(f_origin)) { 
      n0 = length(f_origin); 
      dist_curr = distance(theta_prime = theta_prime, gamma_prime = gamma_prime, n0 = n0, pp=pp, t_unit = t_unit, pad = pad)
      break 
    }
    
    dist_upd = distance(theta_prime = theta_prime, gamma_prime = gamma_prime, n0 = n0, pp=pp, t_unit = t_unit, pad = pad)
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    dist_curr = dist_upd
  }
  
  # print(iter_count)
  
  return(list(n0=n0, dist_min = dist_curr))
}

# Allow n0 to be negative. Move f1 towards f2 by n0. Pad 0 at the head of f1 and f2. 
align_curves_gd = function(f1, f2, n0, step_size, MaxIter=1000, stopping_redu=0.01, pp=FALSE, t_unit=0.05)
{
  f1 = c(rep(0, length(f1)-1), f1) # make N:=length_of_func odd
  f2 = c(rep(0, length(f2)-1), f2)

  
  if (n0 <= 0) {
    r = align_curves_gd_(f_origin = f2, f_shift = f1, n0 = -n0, step_size = step_size, MaxIter = MaxIter, 
                         stopping_redu = stopping_redu, pp=pp, t_unit=t_unit)
    if(r$n0>0)  {return(list(n0 = -r$n0, dist_min = r$dist_min))}
    else {return(align_curves_gd_(f_origin = f1, f_shift = f2, n0 = 0, step_size = step_size, MaxIter = MaxIter, 
                                  stopping_redu = stopping_redu, pp=pp, t_unit=t_unit))}
  }
  else{
    r = align_curves_gd_(f_origin = f1, f_shift = f2, n0 = n0, step_size = step_size, MaxIter = MaxIter, 
                         stopping_redu = stopping_redu, pp=pp, t_unit=t_unit)
    if(r$n0>0) {return(r)}
    else {
      r = (align_curves_gd_(f_origin = f2, f_shift = f1, n0 = 0, step_size = step_size, MaxIter = MaxIter, 
                            stopping_redu = stopping_redu, pp=pp, t_unit=t_unit))
      return(list(n0 = -r$n0, dist_min = r$dist_min))
    }
  }
}


# Initialize n0 by aligning cdf1 and cdf2.
align_pdf_gd = function(pdf1, pdf2, n0=0, step_size=0.02, MaxIter=1000, stopping_redu=0.01, t_unit=0.05)
{
  if(sum(pdf1)==0 || sum(pdf2)==0) n0_init=n0
  else{
    cdf1 = cumsum(pdf1)/sum(pdf1)
    cdf2 = cumsum(pdf2)/sum(pdf2)
    n0_init = align_curves_gd(f1 = cdf1, f2 = cdf2, n0 = n0, step_size = step_size, MaxIter = MaxIter, 
                              stopping_redu = stopping_redu, pp=FALSE, t_unit=t_unit)$n0
  }
  res = align_curves_gd(f1 = pdf1, f2 = pdf2, n0 = n0_init, step_size = step_size, MaxIter = MaxIter, 
                          stopping_redu = stopping_redu, pp=TRUE, t_unit=t_unit)
  return(list(n0 = res$n0, dist_min = res$dist_min))
}



# Align multiple curves with *one* shift parameter --------------------------

get_theta_gamma_prime_v2 = function(f_origin, f_target){
  fft_f_origin = fft(f_origin)
  fft_f_target = fft(f_target)
  N = length(fft_f_origin)
  l_vec = seq(1, N-1)
  pad = tail(f_origin,1)
  
  theta_prime = c(fft_f_origin[1], fft_f_origin[2:N]+1*pad/(1-exp(-1i*2*pi*l_vec/N)))
  gamma_prime = c(fft_f_target[1], fft_f_target[2:N]+1*pad/(1-exp(-1i*2*pi*l_vec/N)))
  
  return(list(theta_prime=theta_prime, gamma_prime=gamma_prime, pad=pad))
}

# Non-negative n0. Move f_origin towards f_shift by n0. If pp=FALSE, pad 1 at the end of f_origin when shifting; otherwise, pad 0.
align_multi_curves_gd_ = function(f_origin_list, f_shift_list, n0, step_size, 
                                  MaxIter=1000, stopping_redu=0.0001, pp=FALSE, t_unit=0.05, weights=NULL)
{
  if(!is.null(weights)&&sum(weights)==0) stop("Weights cannot sum up to zero.")
  if(!is.null(weights)&&sum(weights)!=1) weights = weights/sum(weights)
  
  tmp_list = mapply(FUN = get_theta_gamma_prime, f_origin=f_origin_list, f_shift=f_shift_list, 
                    MoreArgs = list(pp=pp), SIMPLIFY = FALSE)
  theta_prime_list = lapply(X = tmp_list, FUN = "[[", 'theta_prime')
  gamma_prime_list = lapply(X = tmp_list, FUN = "[[", 'gamma_prime')
  pad_list = lapply(X = tmp_list, FUN = "[[", 'pad')
  
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  while (dist_redu>stopping_redu && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd_list = mapply(FUN = gradient, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, pad=pad_list,
                     MoreArgs = list(n0=n0, pp=pp), SIMPLIFY = FALSE)
    if(is.null(weights)) 
      gd = mean(unlist(gd_list))
    else
      gd = sum(unlist(gd_list)*weights)
    
    n0 = n0 - (step_size)*gd
    n0 = round(n0) 
    
    if(n0<0) { 
      n0 = 0; 
      dist_list = mapply(FUN = distance, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, pad=pad_list,
                         MoreArgs = list(n0 = n0, pp=pp, t_unit = t_unit), SIMPLIFY = FALSE)
      if(is.null(weights))
        dist_curr = sqrt(mean(unlist(dist_list)^2))
      else
        dist_curr = sqrt(sum(weights * (unlist(dist_list))^2)) 
      
      break 
    }
    if(n0>length(f_origin_list[[1]])) { 
      n0 = length(f_origin_list[[1]]); 
      # dist_curr = distance(theta_prime = theta_prime, gamma_prime = gamma_prime, n0 = n0, pp=pp, t_unit = t_unit)
      dist_list = mapply(FUN = distance, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, pad=pad_list,
                         MoreArgs = list(n0 = n0, pp=pp, t_unit = t_unit), SIMPLIFY = FALSE)
      if(is.null(weights))
        dist_curr = sqrt(mean(unlist(dist_list)^2))
      else
        dist_curr = sqrt(sum(weights * (unlist(dist_list))^2)) 
      
      break 
    }
    
    # dist_upd = distance(theta_prime = theta_prime, gamma_prime = gamma_prime, n0 = n0, pp=pp, t_unit = t_unit)
    dist_list = mapply(FUN = distance, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, pad=pad_list,
                       MoreArgs = list(n0 = n0, pp=pp, t_unit = t_unit), SIMPLIFY = FALSE)
    if(is.null(weights))
      dist_upd = sqrt(mean(unlist(dist_list)^2))
    else
      dist_upd = sqrt(sum(weights * (unlist(dist_list))^2)) 
    
    
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    dist_curr = dist_upd
  }
  
  # print(iter_count)
  
  return(list(n0=n0, dist_min = dist_curr))
}

# Move f1 towards f2 by n0. Positive n0: towards right. Negative n0: towards left.
align_multi_curves_gd_v2 = function(f_origin_list, f_target_list, n0=0, step_size=0.02, 
                                 MaxIter=1000, stopping_redu=0.0001, t_unit=0.05, weights=NULL,
                                 n0_min = 0, n0_max = length(f_origin_list[[1]]))
{
  if(!is.list(f_origin_list)) f_origin_list = list(f_origin_list)
  if(!is.list(f_target_list)) f_target_list = list(f_target_list)
  
  ### Prevent n0_min and n0_max from further changing
  n0_min -> n0_min
  n0_max -> n0_max 
  
  ### Extend f1 and f2 from [0,T] to [-T,2T]
  extend = function(f) {return(c( rep(head(f,1),length(f)-1), f, rep(tail(f,1),2*round(length(f)/2)) ))} # make N:=length_of_func odd
  f_origin_list = lapply(X = f_origin_list, FUN = extend)
  f_target_list = lapply(X = f_target_list, FUN = extend)
  
  ### Check weights' validity
  if(!is.null(weights)&&sum(weights)==0) stop("Weights cannot sum up to zero.")
  if(!is.null(weights)&&sum(weights)!=1) weights = weights/sum(weights)
  
  
  ### Compute terms needed in gradients
  tmp_list = mapply(FUN = get_theta_gamma_prime_v2, 
                    f_origin=f_origin_list, f_target=f_target_list, 
                    SIMPLIFY = FALSE)
  theta_prime_list = lapply(X = tmp_list, FUN = "[[", 'theta_prime')
  gamma_prime_list = lapply(X = tmp_list, FUN = "[[", 'gamma_prime')
  pad_list = lapply(X = tmp_list, FUN = "[[", 'pad')
  
  ### Gradient descent
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  while (dist_redu>stopping_redu && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd_list = mapply(FUN = gradient_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, pad=pad_list,
                     MoreArgs = list(n0=n0), SIMPLIFY = FALSE)
    if(is.null(weights)) 
      gd = mean(unlist(gd_list))
    else
      gd = sum(unlist(gd_list)*weights)
    
    n0 = n0 - (step_size)*gd
    n0 = round(n0) 
    
    if (n0<n0_min) {
      n0 = n0_min
    }
    if (n0>n0_max) {
      n0 = n0_max
    }
    
    dist_list = mapply(FUN = distance_v2, theta_prime=theta_prime_list, gamma_prime=gamma_prime_list, pad=pad_list,
                       MoreArgs = list(n0 = n0, t_unit = t_unit), SIMPLIFY = FALSE)
    if(is.null(weights))
      dist_upd = sqrt(mean(unlist(dist_list)^2))
    else
      dist_upd = sqrt(sum(weights*(unlist(dist_list))^2)) 
    
    
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    dist_curr = dist_upd
  }
  
  
  if (iter_count==MaxIter) {
    warning("Reached max iteration number when estimating a time shift. Consider adjusting the step size.")
  }
  
  return(list(n0=n0, dist = dist_curr))
  
  
}


### Test
# samp1=c(rnorm(100,5,1), rnorm(30,500,1))
# samp2=c(rnorm(0,25,1), rnorm(200,500,1))
# f1=ecdf(samp1)(seq(0,50,0.1))
# f2=ecdf(samp2)(seq(0,50,0.1))
# plot(f1,type='l',ylim=c(0,1)); lines(f2,col="red")
# res = align_multi_curves_gd_v2(f1,f2,n0_min = 0)
# res$n0
# plot(shift_v2(f1, res$n0),type='l',ylim=c(0,1)); lines(f2,col="red")


