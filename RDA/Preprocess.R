### The default working directory is the current source file location.
### To run the following code, please download the dataset from: 
### https://janelia.figshare.com/s/10833cd5447dbc9aa840
### and save the dataset in '../Zebrafish_spinal_cord_development/'


# Load libraries ----------------------------------------------------------

rm(list=ls())
file_path = "../Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(R.matlab) 
library(data.table)
library(tidyverse)
library(MASS)
library(igraph)
library(ggplot2)
library(grDevices)
library(scales)


# Calculate dF/F values ---------------------------------------------------------------


# We will calculate the dF/F following the method in Wan et al. 
# (Page 25, Step 6)
# Resolution: 4 Hz
# dF/F = (F_cell - F_bsl)/(F_bsl-F_bkgd)
# F_bsl: 20th percentile in the nearby 61 time points
# F_bkgd is set to be zero 

calculate.dFF<-function(F.trace.full,bsl.prob,length.window){
  dFF.full=F.trace.full;
  if(is.nan(F.trace.full[1])){
    
  }else{
    
    for(i in 1:length(F.trace.full)){
      F.trace.ind= max(1,i-(length.window-1)/2):min(length(F.trace.full),i+(length.window-1)/2); 
      F.trace<-F.trace.full[F.trace.ind];
      F.bsl=quantile(F.trace,probs=bsl.prob);  
      dFF.full[i]= (F.trace.full[i]-F.bsl)/F.bsl;  
    }
    
  }
  return(dFF.full)
}

path.list=list.files('../Zebrafish_spinal_cord_development/FunctionalData/');


for(k in 1:length(path.list)){
  path=path.list[[k]]
  dat<- readMat(paste('../Zebrafish_spinal_cord_development/FunctionalData/',path,'/profile.mat',sep=''));
  
  dat.activity=dat$profile.all;
  n.neurons = dim(dat.activity)[1]; 
  n.timeframes=dim(dat.activity)[2];
  
  
  length.windows=61;
  
  # Obtain dFF for each trace
  dat.dFF<-dat.activity;
  for(j in 1:n.neurons){
    dat.dFF[j,]=calculate.dFF(dat.activity[j,],bsl.prob=0.2,length.window=61);
  }
  
  dir.create(path=paste0("../Processed_FunctionalData/", path), recursive = TRUE)
  write.csv(dat.dFF,paste('../Processed_FunctionalData/',path,'/dFF.csv',sep=''),col.names=F)
}






# Basic connectivity inference --------------------------------------------

path.list=list.files('../Processed_FunctionalData/');


rho=0.4; 

window_length = 240 # = 1min at 4 Hz
window_step = 240 # = 1min at 4Hz

# Calculate the covariance matrix 
# Take the inverse of the covariance matrix 

path.list=list.files('../Processed_FunctionalData/');

for(k in 1:length(path.list)){
  path=path.list[[k]]
  dat.dFF=as.matrix(fread(paste('../Processed_FunctionalData/',path,'/dFF.csv',sep='')))
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
  
  adj.full = cor.full.ave;
  for(i in 1:(dim(cor.full.ave)[1])){
    adj.full[i,,] = cor.full.ave[i,,]>rho;
  }
  
  
  ### Find the earliest connecting time for all pairs of neurons
  edge.time = matrix(Inf, nrow=dim(cor.full.ave)[2], ncol=dim(cor.full.ave)[3])
  for(i in 1:nrow(edge.time)){
    for(j in 1:ncol(edge.time)){
      tmp = min(which(adj.full[,i,j]==1)) # index of interval
      edge.time[i,j] = ifelse(tmp<Inf, tmp*(window_length/240),Inf) # edge time (using min as the unit)
    }
  }
  diag(edge.time) = Inf
  
  write.csv(which(!na.inds),paste('../Processed_FunctionalData/',path,'/AvaiNeurons.csv',sep=''))
  
  write.csv(edge.time,paste('../Processed_FunctionalData/',path,'/EdgeTime.csv',sep=''))
  
}




# Extract neural locations and neural biological types --------------------


path.list=list.files('../Zebrafish_spinal_cord_development/FunctionalData/');
for(k in 1:length(path.list)){
  path=path.list[[k]]
  dat=readMat(paste('../Zebrafish_spinal_cord_development/FunctionalData/',path,'/profile.mat',sep=''))
  
  locs.all=cbind(dat$x,dat$y,dat$z)
  mnx = dat$mnx
  
  write.csv(locs.all, paste('../Processed_FunctionalData/',path,'/locs_all.csv',sep=''))
  write.csv(mnx, paste('../Processed_FunctionalData/',path,'/mnx.csv',sep=''))
  
}



