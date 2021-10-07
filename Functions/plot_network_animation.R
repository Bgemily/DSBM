
####### plot network (with igraph)
library(igraph)

plot_network_animation = function(locs, output="./plots/gif/network.gif", 
                                  edge.time=matrix(Inf,nrow=nrow(locs),ncol=nrow(locs)), 
                                  cols=rep(0,nrow(locs)), delay=80, vertex.size=5, alpha=255,
                                  window_list=list(c(0,0)), asp=.5, save_plots=FALSE, remove=TRUE){
  
  T_max = max(edge.time[which(edge.time<Inf)])
  
  # par(mfrow=c(2,1))
  par(mar=c(1,1,1,1))
  
  # get layout of the last network
  window = window_list[[length(window_list)]]
  adj.tmp = edge.time<=min(window[2], T_max) & edge.time>=max(window[1], 2)
  network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected")
  set.seed(16)
  coord = layout_with_fr(network.tmp)
  
  i=100
  time_thres = min(sapply(window_list,"[[",2))
  time_thres = 0
  # time_thres = min(sapply(window_list,"[[",1))
  
  if(save_plots){
    pdf(file=paste0(output,'/', "network_snapshots.pdf"),
        width=1*length(window_list), height=1*asp)
    par(mfrow=c(1,length(window_list)))
    par(mar=c(1.2,0.2,0.2,0.2))
  }
  for (window in window_list) {
    i=i+1
    # if(save_plots){
    #   png(file=paste0(output,'/', "network",i,".png"),
    #       units="in", res=350,
    #       width=2, height=2*asp)
    #   # par(mfrow=c(2,1))
    #   par(mar=c(1.2,0.2,0.2,0.2))
    # }
    for (ax in 2:2 ) {
      adj.tmp = edge.time<=min(window[2], T_max) & edge.time>=max(window[1], 2)
      colnames(adj.tmp) = rownames(adj.tmp) = 1:nrow(edge.time)
      network.tmp = graph_from_adjacency_matrix(adj.tmp, mode="undirected", diag=FALSE)
      
      ends = ends(network.tmp, E(network.tmp))
      if (nrow(ends)>0){
        ends = t(apply(ends, 1, as.numeric))
        edge_cols = apply(ends,1,function(end) ifelse(edge.time[end[1],end[2]]<time_thres,
                                                      rgb(128,128,128,alpha=alpha/10,maxColorValue = 255),
                                                      # rgb(128,128,128,
                                                      #     alpha=alpha/3,maxColorValue = 255) ) )
                                                      rgb(sum(cols[end,1]), sum(cols[end,2]),sum(cols[end,3]),
                                                          alpha=alpha*ifelse(vertex.size[end[1]]==vertex.size[end[2]],
                                                                             yes=1/20, no=1),
                                                          maxColorValue = 255) ) )
        edge.lty = apply(ends,1,function(end) ifelse(end[1]==23 | end[2]==23,
                                                     yes=2, no=1) )
        # edge.lty = apply(ends,1,function(end) 1 )
       
      }
      else{
        edge_cols = NA
        edge.lty = 1
      }
      
      # time_thres = window[2]
      
      
      plot(network.tmp, vertex.label=NA, 
           palette=scales::hue_pal()(3), asp=asp, 
           vertex.size=vertex.size,  
           layout=locs[,c(ax,1)]*cbind(rep(1,nrow(locs)), rep(-1,nrow(locs))),
           edge.color=edge_cols,
           edge.width=0.5,
           edge.lty = edge.lty,
           # edge.color=rgb(0.6,0.8,0.9),
           vertex.frame.color=rgb(cols,maxColorValue=255),
           vertex.color=rgb(cols,maxColorValue=255),
           # vertex.frame.color=rgb(1,0,0),
           # vertex.color=rgb(1,0,0),
           # main=paste(round(window[2]/60),"h",sep='',collapse='~')
      )
      # mtext(paste(round(window[2]/60,digits = 1),"h",sep='',collapse='~'), side=1, line = 0.2)
      mtext(paste(round(window[2],digits = 1),"min",sep='',collapse='~'), side=1, line = 0.2)
      box(col="gray")
    }
    # if(save_plots)
    #   dev.off()
  }
  
  if(save_plots)
    dev.off()
  
  # if(save_plots){
  #   # Converting .png files in one .gif image using ImageMagick
  #   system(paste0("/usr/local/bin/convert -delay ",delay,paste0(" ",output,"/*.png"," ",output,"/"), 'network.gif'))
  # 
  # 
  #   # Remove .png files from working directory
  #   if(remove)
  #     file.remove(list.files(path=output,pattern="*.png",full.names=T))
  # }

  return()
}
