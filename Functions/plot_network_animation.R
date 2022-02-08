
####### Plot network (with igraph)

plot_network_animation = function(locs, node_pair=NA,
                                  vertex.alpha=255,
                                  edge_width_custom=NA,
                                  mar_plot=c(1,0.1,0.1,0.1), box_col="gray",
                                  edge.time=matrix(Inf,nrow=nrow(locs),ncol=nrow(locs)), 
                                  cols=rep(0,nrow(locs)), delay=80, vertex.size=5, alpha=255,
                                  window_list=list(c(0,0)), asp=.5,...){
  
  T_max = max(edge.time[which(edge.time<Inf)])
  
  
  window = window_list[[length(window_list)]]
  adj.tmp = edge.time<=min(window[2], T_max) & edge.time>=max(window[1], 2)
  network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected")
  set.seed(16)
  coord = layout_with_fr(network.tmp)
  
  i=100
  time_thres = min(sapply(window_list,"[[",2))
  p_list = list()
  for (window in window_list) {
    i=i+1
    
    for (ax in 2:2 ) {
      adj.tmp = edge.time<=min(window[2], T_max) & edge.time>=max(window[1], 2)
      colnames(adj.tmp) = rownames(adj.tmp) = 1:nrow(edge.time)
      network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected", diag=FALSE)
      
      ends = ends(network.tmp, E(network.tmp))
      if (nrow(ends)>0){
        ends = t(apply(ends, 1, as.numeric))
        edge_cols = apply(ends,1,function(end) ifelse(edge.time[end[1],end[2]]<time_thres,
                                                      rgb(128,128,128,alpha=alpha/10,maxColorValue = 255),
                                                      rgb(sum(cols[end,1]), sum(cols[end,2]),sum(cols[end,3]),
                                                          alpha=alpha*ifelse(vertex.size[end[1]]==vertex.size[end[2]],
                                                                             yes=1/10, no=1),
                                                          maxColorValue = 255) ) )
        edge.lty = apply(ends,1,function(end) ifelse(!is.na(node_pair) && (node_pair[2] %in% end),
                                                     yes=2, no=1) )
        edge.width = apply(ends,1,function(end) ifelse(!is.na(node_pair) && (node_pair[1] %in% end | node_pair[2] %in% end),
                                                       yes=0.5, no=0.3) )
        
      }
      else{
        edge_cols = NA
        edge.lty = 1
        edge.width = 0.5
      }
      
      if(!is.na(edge_width_custom)){
        edge.width = edge_width_custom
        edge_cols="black"
      }
        
      
      
      par(mar = mar_plot)
      plot(network.tmp, vertex.label=NA, 
           palette=scales::hue_pal()(3), asp=asp, 
           vertex.size=vertex.size,  
           layout=locs[,c(1,ax)],
           edge.color=edge_cols,
           edge.width=edge.width,
           edge.lty = edge.lty,
           vertex.frame.color=rgb(cols,maxColorValue=255),
           vertex.color=rgb(cols, alpha=vertex.alpha, maxColorValue=255),
           ...
      )
      mtext(paste(round(window[2],digits = 1),"min",sep='',collapse='~'), side=1, line = 0.05)
      box(col=box_col)
      
      p = recordPlot()
      p_list = c(p_list, list(p))
    }
    
  }
  
  return(p_list)
}
