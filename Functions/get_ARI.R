# compute ARI for *one* estimated membership
get_one_ARI = function(memb_est_vec, memb_true_vec){
  if (length(memb_est_vec) != length(memb_true_vec)) 
    stop("Two membership vectors should have same length.")
  
  ARI = mclust::adjustedRandIndex(memb_est_vec, memb_true_vec)
  return(ARI)
}


# compute ARI_vec for multiple estimated memberships
get_ARI = function(memb_true_vec, clusters_list=NA, memb_est_vec_list=NA){
  if(is.na(memb_est_vec_list) && is.na(clusters_list))
    stop("At least one of memb_est_vec_list and clusters_list should be non NA.")
  
  if(is.na(memb_est_vec_list))
    memb_est_vec_list = lapply(clusters_list, clus2mem)
  
  ARI_vec = sapply(memb_est_vec_list, get_one_ARI, memb_true_vec)
  return(ARI_vec)
}