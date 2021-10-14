### The default working directory is the current source file location.
### To run the following code, please first download and process the data; see Preprocess.R for instructions.


# Import all functions and packages ----------------------------------------------------

rm(list=ls())
file_path = "../Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(data.table)
library(tidyverse)
library(MASS)
library(igraph)
library(ggplot2)
library(grDevices)
library(scales)



# Sentivitity analysis -----------------------------------------------------------


### Correlation curves (Figure S.3 in SuppMaterial_PartI_DynamicNetworks.pdf) ------

data_folder = "../Processed_FunctionalData/"
path_vec = list.files(data_folder, full.names = TRUE)
path = path_vec[1]

window_length = 240 # = 1min at 4 Hz
window_step = 240 # = 1min at 4Hz


### Read information from data

avai.inds = as.matrix(read.csv(paste(path,'/AvaiNeurons.csv',sep='')))
avai.inds = avai.inds[,-1];

dat.dFF=as.matrix(fread(paste(path,'/dFF.csv',sep='')))
dat.dFF=dat.dFF[,-1]
reduced.dFF=dat.dFF[avai.inds,]


### Calculate the covariance/correlation matrices for each time interval
n.intervals = (dim(reduced.dFF)[2]-window_length) %/% window_step + 1
interval.list = matrix(1:window_length, nrow=n.intervals, ncol=window_length, byrow=TRUE)
interval.list = interval.list + (1:n.intervals-1)*window_step

cor.full=array(dim=c(n.intervals,dim(reduced.dFF)[1],dim(reduced.dFF)[1]))
for(i in 1:n.intervals){
  dat.temp=reduced.dFF[,interval.list[i,]];
  cor.full[i,,]=cor(x=t(dat.temp));  
}

cor.full.ave = cor.full
for(i in 1:(n.intervals-4)){
  cor.full.ave[i,,]=apply(cor.full[i:(min(i+4,n.intervals)),,,drop=F], c(2,3), mean)
}


locs.all = as.matrix(read.csv(paste(path,'/locs_all.csv',sep='')))
locs.all = locs.all[,-1]
locs_mat = locs.all[avai.inds,]



### Visualize correlation curves
id_tmp = which(locs_mat[,2]<0)
plot(1, type="n", xlab="Time (min)", ylab="Correlation", xlim=c(0, 280), ylim=c(-0.1, 1))

for (i in 1:(length(id_tmp)-1)) {
  for (j in (i+1):length(id_tmp)) {
    cor_tmp = cor.full.ave[,id_tmp[i],id_tmp[j]]
    lines(cor_tmp, type='l',col=rgb(0,0,0,0.05),lwd=0.5) 
  }
}

for (i in 2:7) {
  cor_tmp = cor.full.ave[, 1,i]
  lines(cor_tmp, type='l',col=1,lwd=2) 
}








### Type I error with various threshold (Figure S.4 in SuppMaterial_PartI_DynamicNetworks.pdf) -----

set.seed(811)

N_rep = 10

# Simulates a Gaussian process with a given kernel
gaussprocess <- function(from = 0, to = 300, K = function(s, t) {0.8^abs(s-t)*0.03^2/(1-0.8^2)},
                         start = 0, m = 300, mu = rep(0, times = m)) {
  
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = mu, Sigma = Sigma)
  path <- path - path[1] + start  # Must always start at "start"
  
  return(data.frame("t" = t, "xt" = path))
}



# Plots: type I error vs rho, where type I error = (# GP with mean zero & max>rho) / (# GP with mean zero).
rho = 0.8; sigma=0.03
GP_mat = replicate(N_rep, gaussprocess(K=function(s,t)rho^abs(s-t)*sigma^2/(1-rho^2)))
GP_max = apply(GP_mat,2,function(x)max(x$xt))
rho_vec = seq(0.1,0.9,0.1)
type_I_error_vec = rho_vec
for (i in 1:length(rho_vec)) {
  rho = rho_vec[i]
  type_I_error_vec[i] = sum(GP_max>rho) / length(GP_max)
}
plot(y=type_I_error_vec,x=rho_vec,type='b',xlab="Threshold",ylab="Type I error")




# Typical correlation curves and empirical power (Figure S.5 in SuppMaterial_PartI_DynamicNetworks.pdf) ------

set.seed(13)
palette("default")

# Typical correlation curves obtained by applying k-means algorithm
cor_mat = apply(cor.full.ave[,which(locs_mat[,2]<0),which(locs_mat[,2]<0)], 1, function(A)c(A[upper.tri(A)]))
centers = kmeans(cor_mat,5)$centers
order = order(apply(centers,1,max))
centers = centers[order,]

# Visualization of typical correlation curves
plot(1,type='n',xlim=c(0,300),ylim=c(-0.1,1),xlab="Time",ylab="Correlation")
for (k in 1:nrow(centers)) {
  lines(centers[k,],col=k)
}




# Obtain and visualize empirical power under various threshold
plot(1,type='n',xlim=c(min(rho_vec),max(rho_vec)),ylim=c(0,1),xlab="Threshold",ylab="Power")
for (k in 1:nrow(centers)) {
  center = centers[k,]
  GP_mat = replicate(N_rep, gaussprocess(m=length(center),mu=center))
  GP_max = apply(GP_mat,2,function(x)max(x$xt))
  rho_vec = seq(0.1,0.9,0.1)
  power_vec = rho_vec
  for (i in 1:length(rho_vec)) {
    rho = rho_vec[i]
    power_vec[i] = sum(GP_max>rho) / length(GP_max)
  }
  lines(y=power_vec,x=rho_vec,col=k)
}



