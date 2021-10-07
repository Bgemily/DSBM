
library("optparse")

option_list = list(
  make_option(c("-n", "--N_rep"), type="integer", default=5000)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

N_rep = opt$N_rep


# Simulate Gaussian process -------------------------------------------------------------

library(MASS)
gaussprocess <- function(from = 0, to = 300, K = function(s, t) {0.8^abs(s-t)*0.03^2/(1-0.8^2)},
                         start = 0, m = 300, mu = rep(0, times = m)) {
  # Simulates a Gaussian process with a given kernel
  #
  # args:
  #   from: numeric for the starting location of the sequence
  #   to: numeric for the ending location of the sequence
  #   K: a function that corresponds to the kernel (covariance function) of
  #      the process; must give numeric outputs, and if this won't produce a
  #      positive semi-definite matrix, it could fail; default is a Wiener
  #      process
  #   start: numeric for the starting position of the process
  #   m: positive integer for the number of points in the process to simulate
  #
  # return:
  #   A data.frame with variables "t" for the time index and "xt" for the value
  #   of the process
  
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


# explore kernel function ------------------------------------------------------------

# # square exponential kernel
# GP = gaussprocess(from=0,to=300,m=300,K=function(s,t)(0.05)^2*exp(-(s-t)^2/(2*3^2)))
# plot(GP$t,GP$xt, type='l')
# 
# # AR(1) process
# rho = 0.8; sigma=0.03
# GP = gaussprocess(from=0,to=300,m=300,K=function(s,t)rho^abs(s-t)*sigma^2/(1-rho^2))
# plot(GP$t,GP$xt, type='l')
# 
# plot(cor.full.ave[,1,30],type='l')
# 


# type I error ------------------------------------------------------------

set.seed(831)

# N_rep = 10

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
pdf(file="./plots/type_I_error.pdf",width = 4,height = 4)
plot(y=type_I_error_vec,x=rho_vec,type='b',xlab="Threshold",ylab="Type I error")
dev.off()


# white Gaussian kernel
# GP = gaussprocess(from=0,to=300,m=300,K=function(s,t)(0.05)^2*I(s==t))
# plot(GP$t,GP$xt, type='l')



# cluster correlation curves, power, and uncertainty ----------------------------------------------

# N_rep = 10
file = list.files(pattern="cor_full_ave*")[1]
load(file)
file = list.files(pattern="func_20150410*")[1]
load(file)

cor_mat = apply(cor.full.ave[,which(locs[,2]<0),which(locs[,2]<0)], 1, function(A)c(A[upper.tri(A)]))
# centers = cluster::pam(cor_mat,5)$medoids
centers = kmeans(cor_mat,5)$centers
order = order(apply(centers,1,max))
centers = centers[order,]

pdf(file="./plots/typical_correlation.pdf",width = 4,height = 4)
plot(1,type='n',xlim=c(0,300),ylim=c(-0.1,1),xlab="Time",ylab="Correlation")
for (k in 1:nrow(centers)) {
  lines(centers[k,],col=k)
}
dev.off()


# Plots: power vs rho, where power = (# GP with mean=center & max>rho) / (# GP with mean=center).
pdf(file="./plots/power.pdf",width = 4,height = 4)
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
dev.off()

# Plot: uncertainty vs rho, where uncertainty = var(edge time for GP with certain mean function).
pdf(file="./plots/edge_time_uncertainty.pdf",width = 4,height = 4)
plot(1,type='n',xlim=c(min(rho_vec),max(rho_vec)),ylim=c(0,50),xlab="Threshold",ylab="s.d. of edge time")
for (k in 1:nrow(centers)) {
  center = centers[k,]
  GP_mat = replicate(N_rep, gaussprocess(m=length(center),mu=center))
  rho_vec = seq(0.1,0.9,0.1)
  uncertainty_vec = rho_vec
  for (i in 1:length(rho_vec)) {
    rho = rho_vec[i]
    edge_time_vec = apply(GP_mat,2,function(x)which(x$xt>rho)[1])
    uncertainty_vec[i] = sd(edge_time_vec[which(!is.na(edge_time_vec))])
  }
  lines(y=uncertainty_vec,x=rho_vec,col=k)
}
dev.off()
