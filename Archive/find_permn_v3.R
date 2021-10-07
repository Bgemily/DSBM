
### Find the optimal permutation
### Based on v2
### Switch from cdf to fft
### The distance is shift-invariant
### [TO DO: finish the to-do's]
find_permn_v3 = function(center_fft_array_from, center_fft_array_to){
  if (!identical(dim(center_fft_array_from), dim(center_fft_array_to))) {
    stop("dim(center_cdf_array_from) and dim(center_cdf_array_to) should be the same.")
  }
  
  N_clus = dim(center_fft_array_from)[1]
  permn_list = combinat::permn(1:N_clus)
  
  min_dist = Inf
  for (permn in permn_list) {
    
    dist = 0
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        center_fft_1 = center_fft_array_from[permn, permn, ][q,k,]
        center_fft_2 = center_fft_array_to[q,k,]
        
        
        ### [TO DO: update the function to align two Fourier series rather than two cdfs.]
        n0_tmp = align_multi_curves_gd_v2(f_origin_list = list(center_fft_1),
                                          f_target_list = list(center_fft_2),
                                          n0_min = -length(center_fft_1), 
                                          n0_max = length(center_fft_1))$n0
        ### [TO DO: update the shift_v2 function two shift Fourier series rather than cdf/pdf.]
        dist_tmp = sum((center_fft_2-shift_v2(f_origin = center_fft_1,
                                                   n0 = n0_tmp))^2)
        
        
        dist = dist + dist_tmp
      }
    }
  
    ########
    
    if(dist < min_dist){
      the_permn = permn
      min_dist = dist
    }
    
  }
  
  return(list(permn = the_permn, dist = min_dist))
}
