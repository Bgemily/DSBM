
### match clusters so that *clusters are of the same order* across subjects
### F is given as F_true
match_clusters_v3 = function(center_cdf_array_list, true_center_cdf_array){
  N_subj = length(center_cdf_array_list)
  
  ### find permutations by matching clusters of each subject with the first subject
  permn_list = list() 
  for (s in 1:N_subj) {
    res = find_permn_v2(center_cdf_array_from = center_cdf_array_list[[s]], 
                     center_cdf_array_to = true_center_cdf_array)
    permn = res$permn
    permn_list[[s]] = permn
  }
  
  return(list(permn_list=permn_list))
}
