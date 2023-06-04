### The default working directory is the current source file location.

# Load libraries ----------------------------------------------------------

rm(list=ls())
file_path = "./Functions"
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


# Basic connectivity inference --------------------------------------------

subj_name_list = c("func_20150417")
folder_processed_data = '../Processed_FunctionalData/'


rho=0.3; 

window_length = 240 # = 1min at 4 Hz
window_step = 240 # = 1min at 4Hz


for(k in 1:length(subj_name_list)){
  subj_name=subj_name_list[[k]]
  dat.dFF=as.matrix(fread(paste(folder_processed_data,subj_name,'/dFF.csv',sep='')))
  dat.dFF=dat.dFF[,-1]
  # Remove neurons with missing entries
  na.inds= is.na(dat.dFF[,1]);
  reduced.dFF=dat.dFF[!na.inds,];
  
  
  # Calculate the correlation matrices for each time interval
  n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
  interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
  interval.list = interval.list + (1:n.intervals-1)*window_step
  cor.full=array(dim=c(n.intervals,dim(reduced.dFF)[1],dim(reduced.dFF)[1]))
  for(i in 1:n.intervals){
    dat.temp=reduced.dFF[,interval.list[i,]];
    cor.full[i,,]=cor(x=t(dat.temp));  
  }
  
  
  # Threshold the correlation matrix to obtain the connectivity matrix 
  # Nodes are considered as connected if the mean of next 5 correlations is larger than the threshold
  cor.full.ave = cor.full
  for(i in 1:(n.intervals-4)){
    cor.full.ave[i,,]=apply(cor.full[i:(min(i+4,n.intervals)),,,drop=F], c(2,3), mean)
  }
  adj.full = cor.full.ave;
  for(i in 1:(dim(cor.full.ave)[1])){
    adj.full[i,,] = cor.full.ave[i,,]>rho;
  }
  
  
  ### Find the edge time for all pairs of neurons
  edge.time = matrix(Inf, nrow=dim(cor.full.ave)[2], ncol=dim(cor.full.ave)[3])
  for(i in 1:nrow(edge.time)){
    for(j in 1:ncol(edge.time)){
      tmp = min(which(adj.full[,i,j]==1)) # index of interval
      edge.time[i,j] = ifelse(tmp<Inf, tmp*(window_length/240),Inf) # edge time (using min as the unit)
    }
  }
  diag(edge.time) = Inf
  
  write.csv(edge.time,paste(folder_processed_data,subj_name,'/EdgeTime_rho',rho,'.csv',sep=''))
  
}



