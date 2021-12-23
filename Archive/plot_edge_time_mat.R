plot_edge_time_mat = function(edge_time_mat, clusters, n0_mat, denoise=TRUE, 
                              t_vec=seq(0, 50, 0.05), zlim=NULL, reorder=FALSE, showplot=TRUE){
  time_unit = t_vec[2]-t_vec[1]
  
  if(denoise){
    adjs_edge_time_mat = edge_time_mat - time_unit*n0_mat
  }
  else{
    adjs_edge_time_mat = edge_time_mat
  }

  ############### V1
  # diag(adjs_edge_time_mat) = Inf # remove (i,i)
  # 
  # for (q in 1:length(clusters)) {
  #   for (l in 1:length(clusters)) {
  #     
  #     adjs_edge_time_submat = adjs_edge_time_mat[clusters[[q]], clusters[[l]]]
  #     
  #     # NEED JUSTIFICATION
  #     # deal with the situation that shifted event times are negative (thus will be eliminated in estimated cdf)
  #     if (min(adjs_edge_time_submat, na.rm=T) < min(t_vec)) {
  #       MIN_edge_time_submat = min(adjs_edge_time_submat, na.rm=T)
  #       MAX_edge_time_submat = max(adjs_edge_time_submat[is.finite(adjs_edge_time_submat)], na.rm=T)
  #       if (MAX_edge_time_submat-MIN_edge_time_submat <= max(t_vec)-min(t_vec))
  #         adjs_edge_time_submat = adjs_edge_time_submat - MIN_edge_time_submat + min(t_vec)
  #       else{
  #         N_early_events = sum( isTRUE(adjs_edge_time_submat < MAX_edge_time_submat-(max(t_vec)-min(t_vec))) )
  #         N_late_events = sum( (adjs_edge_time_submat > MIN_edge_time_submat+max(t_vec)-min(t_vec)) & is.finite(adjs_edge_time_submat))
  #         if (N_early_events > N_late_events){ # keep lower tail
  #           adjs_edge_time_submat = adjs_edge_time_submat - MIN_edge_time_submat + min(t_vec)
  #           adjs_edge_time_submat[which(adjs_edge_time_submat>max(t_vec))] = Inf
  #         }
  #         else{
  #           adjs_edge_time_submat = adjs_edge_time_submat + max(t_vec) - MAX_edge_time_submat
  #           adjs_edge_time_submat[which(adjs_edge_time_submat<0)] = Inf
  #         }
  #       }
  #     }
  #     
  #     adjs_edge_time_mat[clusters[[q]], clusters[[l]]] = adjs_edge_time_submat
  #     
  #   }
  # }
  # 
  ################ V2
  diag(adjs_edge_time_mat) = Inf # remove (i,i)
  
  # adjs_edge_time_submat = adjs_edge_time_mat
  # # NEED JUSTIFICATION
  # # deal with the situation that shifted event times are negative (thus will be eliminated in estimated cdf)
  # if (min(adjs_edge_time_submat, na.rm=T) < min(t_vec)) {
  #   MIN_edge_time_submat = min(adjs_edge_time_submat, na.rm=T)
  #   MAX_edge_time_submat = max(adjs_edge_time_submat[is.finite(adjs_edge_time_submat)], na.rm=T)
  #   if (MAX_edge_time_submat-MIN_edge_time_submat <= max(t_vec)-min(t_vec))
  #     adjs_edge_time_submat = adjs_edge_time_submat - MIN_edge_time_submat + min(t_vec)
  #   else{
  #     N_early_events = sum( isTRUE(adjs_edge_time_submat < MAX_edge_time_submat-(max(t_vec)-min(t_vec))) )
  #     N_late_events = sum( (adjs_edge_time_submat > MIN_edge_time_submat+max(t_vec)-min(t_vec)) & is.finite(adjs_edge_time_submat))
  #     if (N_early_events > N_late_events){ # keep lower tail
  #       adjs_edge_time_submat = adjs_edge_time_submat - MIN_edge_time_submat + min(t_vec)
  #       adjs_edge_time_submat[which(adjs_edge_time_submat>max(t_vec))] = Inf
  #     }
  #     else{
  #       adjs_edge_time_submat = adjs_edge_time_submat + max(t_vec) - MAX_edge_time_submat
  #       adjs_edge_time_submat[which(adjs_edge_time_submat<0)] = Inf
  #     }
  #   }
  # }
  # adjs_edge_time_mat = adjs_edge_time_submat
  #################################
  
  
  if(reorder){
    adjs_edge_time_mat = adjs_edge_time_mat[unlist(clusters), unlist(clusters)]
  }
  
  if(is.null(zlim)){
    zmax = max(adjs_edge_time_mat[adjs_edge_time_mat<Inf])
    zmin = min(adjs_edge_time_mat)
    zlim = c(zmin,zmax)
  }
  
  if(showplot){
    # fields::image.plot(adjs_edge_time_mat, zlim=zlim, col=fields::tim.colors(300),xaxt="n",yaxt="n")
    colorBreaks = c(seq(1,20,length.out=10),seq(21,100,length.out=300),seq(101,300,length.out=300))
    print(fields::image.plot(adjs_edge_time_mat, zlim=zlim,
                             breaks = colorBreaks,
                             col=fields::tim.colors(length(colorBreaks)-1),
                             # horizontal = TRUE,
                             xaxt="n",yaxt="n"))
  }
  
  # return(adjs_edge_time_mat)
}