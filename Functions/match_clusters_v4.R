
### Match clusters so that *clusters are of the same order* across subjects
### Based on v2
### Switch from cdf to fft
match_clusters_v4 = function(center_fft_array_list){
  N_subj = length(center_fft_array_list)
  
  if (N_subj<=1) {
    stop("Match_clusters() needs at least two subjects.")
  }
  
  ### Find permutations by matching clusters of each subject with the first subject
  permn_list = list() 
  for (s in 1:N_subj) {
    res = find_permn_v3(center_fft_array_from = center_fft_array_list[[s]], 
                        center_fft_array_to = center_fft_array_list[[1]])
    permn = res$permn
    permn_list[[s]] = permn
  }
  
  return(list(permn_list=permn_list))
}
