visual_realdata = function(result, path, suffix="_L.pdf")
{
  tmp = result
  
  # plot estimated f_qk's
  adjs_edge_time_mat = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                                          n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                                          zlim = NULL,reorder = FALSE, showplot = FALSE)
  center_pdf_array = get_center_pdf_array(edge_time_mat = adjs_edge_time_mat, clusters = tmp$clus_result$clusters, 
                                          n0_vec = tmp$clus_result$n0_vec, n0_mat = 0*tmp$clus_result$n0_mat, 
                                          t_vec = tmp$network$t_vec, bw = bw)
  dir.create(paste0("./plots/",path,"/conn_patt/"),recursive = T,showWarnings = F)
  pdf(paste0("./plots/",path,"/conn_patt/",'conn_patt_',path,suffix),width = 4,height = 4)
  g <- plot_pdf_array(center_pdf_array,t_vec = tmp$network$t_vec,y_lim=c(0,0.055))
  print(g)
  dev.off()
  
  # plot edge time matrix heatmap (with or w/out activation time; with or w/out reordering)
  dir.create(paste0("./plots/",path,"/edge_time_mat/"),recursive = T,showWarnings = F)
  pdf(paste0("./plots/",path,"/edge_time_mat/",'edge_time_denoise_reorder_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = TRUE)
  dev.off()
  
  pdf(paste0("./plots/",path,"/edge_time_mat/",'edge_time_noise_reorder_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = TRUE)
  dev.off()
  
  pdf(paste0("./plots/",path,"/edge_time_mat/",'edge_time_denoise_no_order_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = FALSE)
  dev.off()
  
  pdf(paste0("./plots/",path,"/edge_time_mat/",'edge_time_noise_no_order_',path,suffix),width = 4,height = 4)
  plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = FALSE)
  dev.off()
  
  
  ### plot heatmap of matrix of max(vi,vj)
  pdf(paste0("./plots/",path,"/edge_time_mat/",'max_active_time_heatmap_',path,suffix),width = 4,height = 4)
  activ_time_mat = tmp$clus_result$n0_mat * tmp$network$t_vec[2]
  plot_edge_time_mat(edge_time_mat = activ_time_mat, clusters = tmp$clus_result$clusters,
                     n0_mat = tmp$clus_result$n0_mat, denoise = FALSE, t_vec = tmp$network$t_vec, 
                     zlim = c(0,300),reorder = TRUE)
  dev.off()
  
  
  ### plot cluster membership matrix
  dir.create(paste0("./plots/",path,"/clus_mem/"),recursive = T,showWarnings = F)
  pdf(paste0("./plots/",path,"/clus_mem/",'clus_mem_no_order_',path,suffix),width = 4,height = 2)
  image(as.matrix(clus2mem(tmp$clus_result$clusters)), 
        col=gray.colors(length(tmp$clus_result$clusters),start=1,end=0),
        xaxt="n",yaxt='n')
  dev.off()
  
  pdf(paste0("./plots/",path,"/clus_mem/",'clus_mem_reorder_',path,suffix),width = 4,height = 2)
  image(as.matrix(sort(clus2mem(tmp$clus_result$clusters))), 
        col=gray.colors(length(tmp$clus_result$clusters),start=1,end=0),
        xaxt="n",yaxt='n')
  dev.off()
  
  mem_mat = dummies::dummy(clus2mem(tmp$clus_result$clusters))
  pdf(paste0("./plots/",path,"/clus_mem/",'clus_mem_mat_',path,suffix),width = 4,height = 4)
  image(t(mem_mat), col=gray.colors(2,start = 1,end=0), xaxt="n",yaxt='n')
  dev.off()
  
  
  ### plot histogram of activation time
  # clusters_tmp = tmp$clus_result$clusters
  clusters_tmp = tmp$clus_result$clusters[1:(length(tmp$clus_result$clusters)-1)]
  data = data.frame(value=tmp$clus_result$n0_vec[unlist(clusters_tmp)]*tmp$network$t_vec[2], 
                    type=clus2mem(tmp$clus_result$clusters)[unlist(clusters_tmp)])
  dir.create(paste0("./plots/",path,"/active_time/"),recursive = T,showWarnings = F)
  pdf(paste0("./plots/",path,"/active_time/",'active_time_hist_',path,suffix),width = 4,height = 4)
  g<- ggplot(data, aes(x=value,fill=as.factor(type)))+
    geom_histogram( aes(), color="#e9ecef", alpha=0.6, position = 'dodge', bins=5 ) +
    # geom_histogram( aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'dodge', bins=5 ) +
    # geom_density( aes(), color="#e9ecef", alpha=0.6, position = 'dodge', bw=20) +
    scale_fill_manual(values=palette()[2:(length(tmp$clus_result$clusters))]) +
    theme_bw() +
    labs(fill="") #+
    # facet_wrap(~type)
  print(g)
  dev.off()
  
}


visual_realdata_2 = function(result, path, suffix="_L.pdf", edge.time, locs, 
                             reduced.dFF, window_length, window_step, cor.full.ave, member.ship){
  tmp=result
  
  ### plot network development animation (or a snapshot)
  edge.time = tmp$network$edge_time_mat
  edge.time = plot_edge_time_mat(edge_time_mat = tmp$network$edge_time_mat, clusters = tmp$clus_result$clusters,
                                 n0_mat = tmp$clus_result$n0_mat, denoise = TRUE, t_vec = tmp$network$t_vec, 
                                 zlim = NULL,reorder = FALSE, showplot = FALSE)
  
  max_conn_time = max(edge.time[which(edge.time<Inf)])
  mins = lapply(seq(0,max_conn_time,5), function(t)c(0,t))
  
  # plot_network(locs = cbind(locs[tmp$id,2], -locs[tmp$id,1]), edge.time = edge.time, 
  #              output = paste0("animation_denoise",'_',path,suffix, ".gif"),
  #              window_list = list(mins), asp=2, save_plots = T, delay=20, 
  #              cols = t(col2rgb(1+member.ship[tmp$id]*0)), alpha=100)
  
  
  dir.create(paste0("./plots/",path,"/spatial_location/"),recursive = T,showWarnings = F)
  plot_network(locs = locs[,], edge.time = edge.time[,],
               output = paste0('./plots/',path,'/spatial_location/'), 
               filename = paste0("spatial_location_",path,suffix),vertex.size = 3,
               window_list = list(0), asp=0.3, save_plots = T, delay=20,
               cols = t(col2rgb(1+member.ship[])))
  
  ### plot heatmap for neural activity traces
  n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
  interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
  interval.list = interval.list + (1:n.intervals-1)*window_step
  ave_dFF = matrix(nrow=dim(reduced.dFF)[1], ncol=n.intervals)
  for(i in 1:n.intervals){
    ave_dFF[,i] = rowMeans(reduced.dFF[,interval.list[i,]]);
  }
  
  dir.create(paste0("./plots/",path,"/activity_heatmap/"),recursive = T,showWarnings = F)
  pdf(paste0("./plots/",path,"/activity_heatmap/",'activity_heatmap_',path,suffix),width = 4,height = 4)
  fields::image.plot(t(ave_dFF[tmp$id[unlist(tmp$clus_result$clusters)],]),zlim=c(0,0.05), xaxt="n",yaxt='n')
  dev.off()
  
  
  ### plot correlation curves
  dir.create(paste0("./plots/",path,"/corr_curves/"),recursive = T,showWarnings = F)
  pdf(paste0("./plots/",path,"/corr_curves/",'corr_curves_',path,suffix),width = 4,height = 4)
  id_tmp = tmp$id
  plot(1, type="n", xlab="Time (min)", ylab="Correlation", xlim=c(0, 280), ylim=c(-0.1, 1))
  for (i in 1:(length(id_tmp)-1)) {
    for (j in (i+1):length(id_tmp)) {
      cor_tmp = cor.full.ave[,id_tmp[i],id_tmp[j]]
      lines((1:length(cor_tmp))*(window_length/240),cor_tmp, type='l',col=rgb(0,0,0,0.05),lwd=0.5) 
    }
  }
  for (i in 1:6) {
    cor_tmp = cor.full.ave[, sample(id_tmp,1), sample(id_tmp,1)]
    lines((1:length(cor_tmp))*(window_length/240),cor_tmp, type='l',col=i+1,lwd=2) 
  }
  dev.off()
  
  
}










