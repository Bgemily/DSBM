
####### plot network (with igraph)
library(igraph)

plot_network = function(locs, output="./plots/gif/", filename="network.pdf", edge.time=matrix(Inf,nrow=nrow(locs),ncol=nrow(locs)), 
                        cols=rep(0,nrow(locs)), delay=80, vertex.size=5, alpha=255,
                        window_list=list(c(0,0)), asp=.5, save_plots=FALSE, remove=TRUE){
  
  T_max = max(edge.time[which(edge.time<Inf)])
  
  .pardefault = par(no.readonly = T)
  # par(mfrow=c(2,1))
  par(mar=c(1,1,1,1))
  
  # get layout of the last network
  window = window_list[[length(window_list)]]
  adj.tmp = edge.time<=min(window[2], T_max) & edge.time>=max(window[1], 2)
  network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected")
  # set.seed(16)
  coord = layout_with_fr(network.tmp)
  
  i=100
  for (window in window_list) {
    i=i+1
    
    for (ax in 2:3 ) {
      if(save_plots){
        pdf(file=paste0(output,'ax_',ax,"_",filename), 
            width=3, height=3*asp)
        # par(mfrow=c(2,1))
        par(mar=c(1,1,1,1))
      }
      adj.tmp = edge.time<=min(window[2], T_max) & edge.time>=max(window[1], 2)
      colnames(adj.tmp) = rownames(adj.tmp) = 1:nrow(edge.time)
      network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected", diag=FALSE)
      
      ends = ends(network.tmp, E(network.tmp))
      if (nrow(ends)>0){
        ends = t(apply(ends, 1, as.numeric))
        edge_cols = apply(ends,1,function(end) rgb(mean(cols[end,1]), mean(cols[end,2]),mean(cols[end,3]), 
                                                   alpha=alpha,maxColorValue = 255))
      }
      else{
        edge_cols = NA
      }
        
      
      plot(network.tmp, vertex.label=NA, 
           palette=scales::hue_pal()(3), asp=asp, 
           vertex.size=vertex.size,  
           layout=locs[,c(1,ax)],
           edge.color=edge_cols,
           # edge.color=rgb(0.6,0.8,0.9),
           vertex.frame.color=rgb(cols,maxColorValue=255),
           vertex.color=rgb(cols,maxColorValue=255),
           # vertex.frame.color=rgb(1,0,0),
           # vertex.color=rgb(1,0,0), 
           # main=paste(round(window[2]),"min",sep='',collapse='~')
           )
      # title(paste(round(window[2]),"min",sep='',collapse='~'), line = -1, adj=1)
      box(col="gray")
      
      if(save_plots)
        dev.off()
      
    }
    
  }
  
  
  par(.pardefault)
  return()
}
