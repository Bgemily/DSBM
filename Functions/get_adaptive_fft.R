### Choose adaptive frequency truncation parameter to smooth a vector of event times
get_adaptive_fft = function(event_time_vec, freq_trun_max, t_vec, freq_trun_min=freq_trun_max){
  if (length(t_vec)<2*freq_trun_max+1) {
    stop("Length of t_vec should be no less than 2*freq_trun_max+1.")
  }
  
  loss_min = Inf
  ### Get empirical intensity of event times
  emp_intens_vec = hist(event_time_vec, breaks=t_vec, plot=FALSE)$counts / length(event_time_vec)
  emp_intens_vec = c(0,emp_intens_vec)
  ### Get normalized fourier series
  emp_intens_fft = fft(emp_intens_vec) / length(t_vec)
  ### Get truncated fourier series
  fft_vec_max = c(tail(emp_intens_fft, freq_trun_max),
                  head(emp_intens_fft, freq_trun_max+1))
  
  for (freq_trun in freq_trun_min:freq_trun_max) {
    ### Get fft_vec
    fft_vec = c(tail(emp_intens_fft, freq_trun),
                head(emp_intens_fft, freq_trun+1))
    
    # Calculate loss (c.f. Bigot and Jendre'13, Eq.(11))
    fft_vec_extend = c(rep(0,freq_trun_max-freq_trun),
                       fft_vec,
                       rep(0,freq_trun_max-freq_trun))
    bias = sum(abs(fft_vec_max-fft_vec_extend)^2)
    variance = 2*(2*freq_trun+1)*(sum(event_time_vec<max(t_vec))/(length(event_time_vec))^2)
    loss = bias + variance
    
    # Update best frequency truncation and fft
    if (loss < loss_min) {
      loss_min = loss
      freq_trun_best = freq_trun
      fft_vec_best = fft_vec_extend
    }
  }
  
  return(list(freq_trun_best = freq_trun_best, 
         fft_vec_best = fft_vec_best))
}