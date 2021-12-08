
##### plot local traces

library(ggridges)
library(ggplot2)
library(reshape2)
library(gridExtra)

plot_local_traces = function(locs, reduced.dFF, edge_time_mat=NULL, snapshot_mins=c(30,120,240), 
                             membership=1,
                             window_length=1, id_vec=NULL, show_plot=FALSE,  scale=2){
  order.temp = order(locs[,1], decreasing = TRUE)
  order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])
  
  order = c(which(membership[order.temp]==0),
            which(membership[order.temp]==3),
            which(membership[order.temp]==2),
            which(membership[order.temp]==1))
  membership = membership[order.temp][order]
  # reduced.dFF = reduced.dFF[order.temp,][order,]
  
  
  if(!is.null(edge_time_mat)){
    tmp_ind = 1:nrow(locs)
    min_edge_time = sapply(1:length(tmp_ind), function(i)min(edge_time_mat[tmp_ind,tmp_ind][i,-i]))
    order.temp = order(min_edge_time); 
    order.temp = rev(order.temp)
    order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])
  }
  
  
  
  if(!is.null(id_vec))
    id_vec = sapply(id_vec, function(id)which(order.temp==id))
  
  timefield = lapply(snapshot_mins, function(min)min*240+1:(240*window_length))
  g.temp = list()
  for (t in timefield) {
    mat.temp = reduced.dFF[,t]

    df.temp = data.frame(id=1:nrow(mat.temp), membership=membership, mat.temp)
    df.temp = pivot_longer(df.temp, cols=!c(id, membership), names_to = "time", values_to = "dFF")
    g = df.temp %>% 
      mutate(membership=as.factor(membership)) %>%
      ggplot(aes(x = time, y = id, height = dFF, group=id, color=membership),fill="white") +
      geom_ridgeline(scale=scale, min_height=-0.1, alpha=0, size=0.5) +
      scale_color_manual(values = palette()[c(8,2:4)], aesthetics = c('color'))+
      xlab(paste(round((t[1]-1)/240,digits = 1), "min"))+
      # theme_bw()+
      theme(legend.position = "none", axis.text.x=element_blank(), 
            axis.ticks.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank(), axis.title.y = element_blank()) 
    g.temp = c(g.temp, list(g))
  }
  if(show_plot)
    do.call("grid.arrange", c(g.temp, ncol=length(g.temp)))
  else
    return(g.temp)
}
