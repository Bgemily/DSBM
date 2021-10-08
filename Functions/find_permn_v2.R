
### Find the optimal permutation
### The distance is shift-invariant
### THe distance considers both probability and shape
find_permn_v2 = function(center_cdf_array_from, center_cdf_array_to){
  if (!identical(dim(center_cdf_array_from), dim(center_cdf_array_to))) {
    stop("dim(center_cdf_array_from) and dim(center_cdf_array_to) should be the same.")
  }
  
  N_clus = dim(center_cdf_array_from)[1]
  permn_list = combinat::permn(1:N_clus)
  
  min_dist = Inf
  for (permn in permn_list) {
    
    dist = 0
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        center_cdf_1 = center_cdf_array_from[permn, permn, ][q,k,]
        center_cdf_2 = center_cdf_array_to[q,k,]
        if (var(center_cdf_1)==0) {
          center_cdf_1 = c(center_cdf_1[-1],center_cdf_1[1]+1e-10)
        }
        if (var(center_cdf_2)==0) {
          center_cdf_2 = c(center_cdf_2[-1],center_cdf_2[1]+1e-10)
        }
        # tmp = 1-max(ccf(x = center_cdf_1, y = center_cdf_2, plot=FALSE)$acf)
        
        conn_prob_1 = tail(center_cdf_1, 1)
        conn_prob_2 = tail(center_cdf_2, 1)
        center_cdf_1_normed = center_cdf_1 / conn_prob_1
        center_cdf_2_normed = center_cdf_2 / conn_prob_2
        
        n0_tmp = align_multi_curves_gd_v2(f_origin_list = list(center_cdf_1_normed),
                                          f_target_list = list(center_cdf_2_normed),
                                          n0_min = -length(center_cdf_1_normed), 
                                          n0_max = length(center_cdf_1_normed))$n0
        dist_1 = sum((center_cdf_2_normed-shift_v2(f_origin = center_cdf_1_normed,
                                                   n0 = n0_tmp))^2)
        dist_2 = (conn_prob_2-conn_prob_1)^2
        dist_tmp = dist_1 * (conn_prob_1^2+conn_prob_2^2)/2 + 
                    dist_2 * ( sum((center_cdf_2_normed*I(center_cdf_2_normed<1))^2)+
                               sum((center_cdf_1_normed*I(center_cdf_1_normed<1))^2) )/2
        
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
