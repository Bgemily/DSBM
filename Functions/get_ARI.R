# Compute ARI for estimated membership
get_one_ARI = function(memb_est_vec, memb_true_vec){
  if (length(memb_est_vec) != length(memb_true_vec)) 
    stop("Two membership vectors should have same length.")
  
  ARI = mclust::adjustedRandIndex(memb_est_vec, memb_true_vec)
  return(ARI)
}

