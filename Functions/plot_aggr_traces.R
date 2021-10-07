plot_aggr_traces = function(edge_time_mat, xlim=c(0,50), ylim=c(0,0.1), bw=NULL, color=NULL){
  i=1; plot(0, lwd=0, xlim=xlim,ylim=ylim, xlab="time(min)", main="Event Rate", ylab=NA)
  for(i in 1:nrow(edge_time_mat)){
    if(length(which(edge_time_mat[i,-i]<Inf))>=2){
      if(!is.null(bw))
        density = density(edge_time_mat[i,-i], bw=bw, from=xlim[1], to=xlim[2])
      else
        density = density(edge_time_mat[i,-i], from=xlim[1], to=xlim[2])
      
      col = ifelse(is.null(color),1,color[i])
      lines(density$x, density$y * (sum(edge_time_mat[i,-i]<Inf))/(nrow(edge_time_mat)-1),col=col)
    }
  }
  
}