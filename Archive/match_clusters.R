
### match clusters so that *clusters are of the same order* across subjects
match_clusters = function(center_cdf_array_list){
  N_subj = length(center_cdf_array_list)
  
  if (N_subj<=1) {
    stop("Match_clusters() needs at least two subjects.")
  }
  
  ### find permutations by matching clusters of each subject with the first subject
  permn_list = list() 
  for (s in 1:N_subj) {
    res = find_permn(center_cdf_array_from = center_cdf_array_list[[s]], 
                     center_cdf_array_to = center_cdf_array_list[[1]])
    permn = res$permn
    permn_list[[s]] = permn
  }
  
  return(list(permn_list=permn_list))
}
