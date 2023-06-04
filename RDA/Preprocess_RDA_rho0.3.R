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
# Input: folder_processed_data, subj_name_list, rho
# Output: edge.time

folder_processed_data = '../Processed_FunctionalData/'
subj_name_list = c("func_20150417")
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


# Save index of neurons to be analyzed in left / right spines (based on edge times obtained from rho=0.4) ------------------------------------
# Input: folder_processed_data, subj_name_list
# Output: inds_neuron_analyzed_L, inds_neuron_analyzed_R

folder_processed_data = "../Processed_FunctionalData/"
subj_name_list = c("func_20150417")
avai_inds_list = list()
for(m in 1:length(subj_name_list)){ 
  subj_name = subj_name_list[m]
  
  ### Read information from data
  edge_time_mat = as.matrix(read.csv(paste(folder_processed_data,subj_name, '/EdgeTime.csv', sep='')))
  edge_time_mat = edge_time_mat[,-1]
  avai.inds = as.matrix(read.csv(paste(folder_processed_data,subj_name,'/AvaiNeurons.csv',sep='')))
  avai.inds = avai.inds[,-1];
  locs.all = as.matrix(read.csv(paste(folder_processed_data,subj_name,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs_mat = locs.all[avai.inds,]
  ### Separate networks in left and right spines
  inds_L = which(locs_mat[,2]<0)
  inds_R = which(locs_mat[,2]>0)
  locs_mat_L = locs_mat[inds_L,]
  locs_mat_R = locs_mat[inds_R,]
  edge_time_mat_L = edge_time_mat[inds_L, inds_L]
  edge_time_mat_R = edge_time_mat[inds_R, inds_R]
  ### Remove neurons with no connection and neurons with abnormal behavior
  avai_inds_L = which(rowSums(edge_time_mat_L<Inf)>0 )
  avai_inds_R = which(rowSums(edge_time_mat_R<Inf)>0 )
  if (TRUE | m==2) {
    avai_inds_L = avai_inds_L[-c(1,2)]
    avai_inds_R = avai_inds_R[-c(45,34,15)]
  }
  inds_neuron_analyzed_L = inds_L[avai_inds_L]
  inds_neuron_analyzed_R = inds_R[avai_inds_R]

  write.csv(inds_neuron_analyzed_L, paste(folder_processed_data,subj_name,'/Index_neuron_analyzed_L.csv',sep=''))
  write.csv(inds_neuron_analyzed_R, paste(folder_processed_data,subj_name,'/Index_neuron_analyzed_R.csv',sep=''))
  
}

