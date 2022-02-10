est_edge_time = function(window_length=240, window_step=window_length, rho=0.6, 
                         processed_data_folder='../processed_FunctionalData/', path){
  
  library(data.table)


  dat.dFF=as.matrix(fread(paste(processed_data_folder,path,'/dFF.csv',sep='')))
  dat.dFF=dat.dFF[,-1]
  # Simply calculate the covariance among dFF within the time periods
  # Obtain the list of missing traces 
  na.inds= is.na(dat.dFF[,1]);
  reduced.dFF=dat.dFF[!na.inds,];
  
  # Calculate the covariance/correlation matrices for each time interval
  
  
  n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
  interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
  interval.list = interval.list + (1:n.intervals-1)*window_step
  
  
  cor.full=array(dim=c(n.intervals,dim(reduced.dFF)[1],dim(reduced.dFF)[1]))
  for(i in 1:n.intervals){
    dat.temp=reduced.dFF[,interval.list[i,]];
    cor.full[i,,]=cor(x=t(dat.temp));  
  }
  
  
  
  # Thresholding the correlation matrix to obtain the connectivity matrix 
  # nodes are considered as connected if the mean of next 5 correlations is larger than the threshold
  
  cor.full.ave = cor.full
  for(i in 1:(n.intervals-4)){
    cor.full.ave[i,,]=apply(cor.full[i:(min(i+4,n.intervals)),,,drop=F], c(2,3), mean)
  }
  adj.full=cor.full;
  for(i in 1:(n.intervals-4)){
    adj.full[i,,]=cor.full.ave[i,,]>rho;
  }
  
  
  # Find the minimum connecting time for the neurons 
  edge.time=matrix(Inf,nrow=dim(reduced.dFF)[1],ncol=dim(reduced.dFF)[1])
  for(i in 1:dim(reduced.dFF)[1]){
    for(j in 1:dim(reduced.dFF)[1]){
      tmp = min(which(adj.full[,i,j]==1)) # index of interval
      edge.time[i,j] = ifelse(tmp<Inf, interval.list[tmp,1]/240,Inf) # edge time (using min as the unit)
    }
  }
  
  
  write.csv(edge.time,paste(processed_data_folder,path,'/EdgeTime.csv',sep=''),col.names=F)
  save(list = c("cor.full.ave"),
       file=paste0(processed_data_folder,path,'/cor_full_ave','_win',window_length,'_rho',rho,'.rdata'))

  return(list(edge_time_mat = edge.time, cor.full.ave = cor.full.ave))
  
}



