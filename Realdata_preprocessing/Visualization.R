
setwd("/Users/bgemily/Documents/Academic/SC/graphon/simulation")

library(data.table)
library(R.matlab)


path.list=list.files('../processed_FunctionalData/');



# Plot neuron activities at certain time points --------------------------------


library(ggplot2)
for (k in 1:length(path.list)) {
  path=path.list[[k]]

  edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
  edge.time=edge.time[,-1]
  
  avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
  avai.inds=avai.inds[,-1];
  
  # member.ship = as.matrix(read.csv(paste('./processed_FunctionalData/',path,'/MembShip.csv',sep='')))
  # member.ship=member.ship[,-1];
  
  dat.dFF=as.matrix(fread(paste('../processed_FunctionalData/',path,'/dFF.csv',sep='')))
  dat.dFF=dat.dFF[,-1]
  reduced.dFF=dat.dFF[avai.inds,];
  
  locs.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs=locs.all[avai.inds,]
  
  mnx.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/mnx.csv',sep='')))
  mnx.all = mnx.all[,-1]
  mnx = mnx.all[avai.inds]
  
  # islet = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/islet.csv',sep='')))
  # islet = islet[,-1]
  
  N_timepoints = ncol(reduced.dFF)
  mins = c(30,120,240)
  if (N_timepoints<57600) # 4h = 57600 time points
    mins = c(30, 60, 120, round(N_timepoints/240-5))
  # g.temp = plot_local_traces(locs=locs, reduced.dFF=reduced.dFF, snapshot_mins=mins, window_length=1)
  g.temp = plot_local_traces(locs=locs[tmp$id,], reduced.dFF=reduced.dFF[tmp$id,], snapshot_mins=mins, window_length=1,scale=4)
  g.temp = gridExtra::arrangeGrob(grobs = g.temp, ncol = length(g.temp))
  grid::grid.draw(g.temp)
  ggsave(filename = paste0('./plots/',path,'_local_traces.pdf'), plot = g.temp,
         width = 6, height = 3, units = "in")
}



# Plot neuron activities for whole development process -------------------------


for (t in seq(282,0,-3)) {
  g_tmp = plot_local_traces(locs[tmp_ind,],reduced.dFF[tmp_ind,], edge_time_mat[tmp_ind,tmp_ind],  snapshot_mins = c(t),window_length = 3)
  g_tmp = gridExtra::arrangeGrob(grobs = g_tmp, ncol = length(g_tmp))
  ggsave(filename = paste0('./plots/',path,'_traces_','Rside_',t,'.pdf'), plot = g_tmp,
         width = 8.4, height = 5.35, units = "in")
}



# plot network development process ----------------------------------------


for (k in 1:length(path.list)) {
  path=path.list[[k]]
  
  edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
  edge.time=edge.time[,-1]
  
  avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
  avai.inds=avai.inds[,-1];
  
  member.ship = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/MembShip.csv',sep='')))
  member.ship=member.ship[,-1];
  
  locs.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs=locs.all[avai.inds,]
  
  mnx.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/mnx.csv',sep='')))
  mnx.all = mnx.all[,-1]
  mnx = mnx.all[avai.inds]
  
  # mins = list(c(0,120), c(0,150), c(0,180))
  max_conn_time = max(edge.time[which(edge.time<Inf)])
  mins = lapply(seq(60,max_conn_time+15,5), function(t)c(0,t))
  
  L_id = which(locs[,2]<0 & rowSums(edge.time[,which(locs[,2]<0)]<Inf)>=2)
  R_id = which(locs[,2]>0 & rowSums(edge.time[,which(locs[,2]>0)]<Inf)>=2)
  tmp_id = L_id
  # plot_network(locs = locs[tmp_id,], edge.time = edge.time[tmp_id,tmp_id], 
  #              output = "Lside.gif",
  #              window_list = mins, asp=1, save_plots = T, delay=20)
  plot_network(locs = cbind(locs[tmp_id,2], -locs[tmp_id,1]), edge.time = edge.time[tmp_id,tmp_id], 
               output = "animation.gif",
               window_list = mins, asp=2, alpha=100,
               save_plots = T, delay=20, 
               cols = t(col2rgb(1+member.ship[tmp_id]*0)))
  
  plot_network(locs = locs[,], edge.time = edge.time[,], 
               output = "Network.pdf", vertex.size = 3,
               window_list = list(0), asp=0.3, save_plots = T, delay=20, 
               cols = t(col2rgb(1+member.ship[])))
  
  
  order.temp = order(locs[,1], decreasing = TRUE)
  order.temp = c(order.temp[locs[order.temp,2]<0], order.temp[locs[order.temp,2]>0])
  
}



# plot edge time matrix heatmap -------------------------------------------

for (k in 1:length(path.list)) {
  path=path.list[[k]]
  
  edge.time=as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/EdgeTime.csv',sep='')))
  edge.time=edge.time[,-1]
  
  avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
  avai.inds=avai.inds[,-1];
  
  member.ship = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/MembShip.csv',sep='')))
  member.ship=member.ship[,-1];
  
  locs.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs=locs.all[avai.inds,]
  
  mnx.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/mnx.csv',sep='')))
  mnx.all = mnx.all[,-1]
  mnx = mnx.all[avai.inds]
  
  
  L_id = which(locs[,2]<0 & rowSums(edge.time[,which(locs[,2]<0)]<Inf)>=2)
  R_id = which(locs[,2]>0 & rowSums(edge.time[,which(locs[,2]>0)]<Inf)>=2)
  tmp_id = R_id
  tmp_id = tmp_id[order(member.ship[tmp_id])]

  fields::image.plot(edge.time[tmp_id,tmp_id],zlim=c(0,300), col=fields::tim.colors())

}

load(list.files(pattern=path)) #load .rdata
tmp = L_result
edge_time_mat = tmp$network$edge_time_mat
diag(edge_time_mat) = Inf
tau_mat = tmp$clus_result$n0_mat*tmp$network$t_vec[2]
edge_time_mat_denoise = edge_time_mat - tau_mat
diag(edge_time_mat_denoise) = Inf
edge_time_mat_denoise = edge_time_mat_denoise - min(edge_time_mat_denoise)





#  Plot edge time's density  ----------------------------------------------------------------------


image(edge.time)

bw=10
edge.time.tmp = edge_time_mat
i=2; plot(density(edge.time.tmp[i,-i])$x, density(edge.time.tmp[i,-i])$y * (sum(edge.time.tmp[i,-i]<Inf))/(nrow(edge.time.tmp)-1), 
          lwd=0, type='l',xlim=c(0,340),ylim=c(0,0.004), xlab="time(min)", main="Event Rate", ylab=NA)
for(i in tmp_ind[1:30]){
  if(length(which(edge.time.tmp[i,-i]<Inf))>=2)
    density = density(edge.time.tmp[i,-i], bw=bw)
    lines(density$x, density$y * (sum(edge.time.tmp[i,-i]<Inf))/(nrow(edge.time.tmp)-1),
          xlim=c(0,340),col=i)
}


edge.time.tmp = edge.time
i=1; plot(ecdf(edge.time.tmp[i,-i])(1:340), 
          lwd=0, type='l',xlim=c(0,340),ylim=c(0,0.6), xlab="time(min)", main="Event Rate", ylab=NA)
for(i in 1:nrow(edge.time.tmp)){
  if(length(which(edge.time.tmp[i,-i]<Inf))>=2)
    density = density(edge.time.tmp[i,-i], bw=10)
  lines(ecdf(edge.time.tmp[i,-i])(1:340), 
        xlim=c(0,340), lwd=0.3)
}



# Animation of neuron movement  --------------------------------------------------------------

library(gganimate)
path = path.list[[7]]
track_smooth_all = readRDS(paste('../processed_FunctionalData/',path,'/tracks_smoothed.rds',sep=''))
avai.inds = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/AvaiNeurons.csv',sep='')))
avai.inds=avai.inds[,-1];

mnx.all = as.matrix(read.csv(paste('../processed_FunctionalData/',path,'/mnx.csv',sep='')))
mnx.all = mnx.all[,-1]
mnx = mnx.all[avai.inds]


track_smooth = track_smooth_all[avai.inds, seq(1,dim(track_smooth_all)[2],1200),]

N_node = dim(track_smooth)[1]
N_timepoints = dim(track_smooth)[2]

track_smooth_df = matrix(track_smooth, N_node*N_timepoints, dim(track_smooth)[3])
track_smooth_df = data.frame(x=track_smooth_df[,1],y=track_smooth_df[,2], z=track_smooth_df[,3],
           time=rep(1:N_timepoints,each=N_node), id=rep(1:N_node,time=N_timepoints))

anim = ggplot(track_smooth_df, aes(x=x,y=y,col=mnx[id]))+
        geom_point()+
        transition_time(time) 
animate(anim, renderer = gifski_renderer())
# anim_save("func_20150410_xy.gif")

anim = ggplot(track_smooth_df, aes(x=x,y=z,col=mnx[id]))+
  geom_point()+
  transition_time(time) 
animate(anim, renderer = gifski_renderer())
# anim_save("func_20150410_xz.gif")




# explore the data --------------------------------------------------------

edge_time_mat[tmp_ind,tmp_ind][1:20,1:20]

plot(edge_time_mat[tmp_ind,tmp_ind][1,],type='b',cex=.3,ylim=c(0,300))
lines(edge_time_mat[tmp_ind,tmp_ind][3,],type='b',cex=.3,col=2)
lines(edge_time_mat[tmp_ind,tmp_ind][4,],type='b',cex=.3,col=3)
lines(edge_time_mat[tmp_ind,tmp_ind][33,],type='b',cex=.3,col=4)
lines(edge_time_mat[tmp_ind,tmp_ind][34,],type='b',cex=.3,col=5)
lines(edge_time_mat[tmp_ind,tmp_ind][35,],type='b',cex=.3,col=6)
lines(edge_time_mat[tmp_ind,tmp_ind][25,],type='b',cex=.3,col=7)
lines(edge_time_mat[tmp_ind,tmp_ind][26,],type='b',cex=.3,col=8)




i=11;
plot(cor.full[,tmp_ind[13],tmp_ind[i]],type='l',col=1,xlim=c(0,250),ylim=c(0,1)); 
lines(cor.full[,tmp_ind[9],tmp_ind[i]],type='l',col=2)
lines(cor.full[,tmp_ind[22],tmp_ind[i]],col="green")
lines(cor.full[,tmp_ind[24],tmp_ind[i]],col=4)
lines(cor.full[,tmp_ind[25],tmp_ind[i]],col=5)

plot(cor.full[,23,2],type='l'); 
col=1
for (i in c(4,6,8)) {
  lines(cor.full[,34,tmp_ind[i]],col=rainbow(5)[col])
  col=col+1
}
lines(cor.full[,tmp_ind[4],tmp_ind[8]],col='red')

t=111
i=c(10,12)
g_tmp = plot_local_traces(locs[rep(1,length(tmp_ind[i])),],reduced.dFF[tmp_ind[i],], edge_time_mat[tmp_ind[i],tmp_ind[i]],  snapshot_mins = c(t),window_length = 3)
g_tmp
cov(reduced.dFF[tmp_ind[1],t*240+1:(1*240)],reduced.dFF[tmp_ind[2],t*240+1:(1*240)])


magnitude = c()
for (t in seq(0,280,1)) {
  tmp = apply(reduced.dFF[,t*240+1:(1*240)], 1, max)
  magnitude = cbind(magnitude,tmp)
}

for (i in L_result$id[L_result$clus_result$clusters[[3]]]) {
  plot(magnitude[i,],type='l',main=i,col=2)
}
plot(magnitude[16,],type='l')
lines(magnitude[4,],type='l',col="blue")


t=50
magnitude = apply(reduced.dFF[tmp_ind,t*240+1:(3*240)], 1, max)
which(magnitude>0.2)
plot(log(magnitude),(rowMeans(cor.full[t,tmp_ind,tmp_ind])))



