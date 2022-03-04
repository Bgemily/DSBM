### The default working directory is the current source file location.
### To run the following code, please first download and process the data; see Preprocess.R for instructions.


# Import all functions and packages ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(data.table)
library(tidyverse)
library(MASS)
library(igraph)
library(ggplot2)
library(grDevices)
library(scales)
library(forecast)



# Sentivitity analysis -----------------------------------------------------------


### Correlation curves (Figure S.3 in SuppMaterial_PartI_DynamicNetworks.pdf) ------

data_folder = "../Processed_FunctionalData/"
path_vec = list.files(data_folder, full.names = TRUE)
path = path_vec[2]

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



### Visualize correlation curves ####
##### Left ####
data_1_left = data.frame()
id_tmp = which(locs_mat[,2]<0)
cor.full.ave_left = cor.full.ave[,id_tmp,id_tmp]
for (i in 1:(length(id_tmp)-1)) {
  for (j in (i+1):length(id_tmp)) {
    cor_tmp = cor.full.ave_left[ ,i,j]
    tmp = bind_cols(correlation=cor_tmp,
                    time=1:length(cor_tmp),
                    node_1=i,
                    node_2=j)
    data_1_left = bind_rows(data_1_left, tmp)
  }
}

set.seed(1329)
data_1_left = data_1_left %>%
  group_by(node_1,node_2) %>%
  mutate(highlight = rep(sample(0:1,1,prob=c(0.996,0.004))*(max(correlation)>0.2), n()))
data_1_left[which(data_1_left$node_1==2 & data_1_left$node_2==46), 'highlight'] = 2
g1_left = data_1_left %>%  
  mutate(highlight = factor(highlight, levels = c(0,1,2))) %>%
  ggplot() +
  geom_line(mapping = aes(x=time, y=correlation, 
                          group=interaction(node_1,node_2,highlight),
                          color=highlight,
                          alpha=highlight,
                          size=highlight),
            ) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('darkgray','black','red'),
                     breaks = c(0,1,2)) +
  scale_alpha_manual(values=c(0.07,0.7,0.7),
                     breaks = c(0,1,2)) +
  scale_size_manual(values=c(0.2, 0.5,0.5),
                    breaks = c(0,1,2)) +
  xlab("Time (min)") +
  ylab("Correlation") +
  scale_y_continuous( breaks = seq(0,1,0.2)) +
  coord_cartesian(xlim=c(0, 330), ylim=c(-0.1, 1))
# g1_left

pdf(file=paste0("../Results/Plots/Temp/Sensitivity/",
                "Correlation_curves_L",
                ".pdf"),
    width = 6, height = 2)
print(g1_left)
dev.off()

##### Right ####
data_1_right = data.frame()
id_tmp = which(locs_mat[,2]>0)
cor.full.ave_right = cor.full.ave[,id_tmp,id_tmp]
for (i in 1:(length(id_tmp)-1)) {
  for (j in (i+1):length(id_tmp)) {
    cor_tmp = cor.full.ave_right[ ,i,j]
    tmp = bind_cols(correlation=cor_tmp,
                    time=1:length(cor_tmp),
                    node_1=i,
                    node_2=j)
    data_1_right = bind_rows(data_1_right, tmp)
  }
}

set.seed(8119)
data_1_right = data_1_right %>%
  group_by(node_1,node_2) %>%
  mutate(highlight = rep(sample(0:1,1,prob=c(0.997,0.003))*(max(correlation)>0.2), n())) %>%
  mutate(highlight = factor(highlight, levels=c(0,1)))
g1_right = data_1_right %>%
  ggplot() +
  geom_line(mapping = aes(x=time, y=correlation, 
                          group=interaction(node_1,node_2,highlight),
                          color=highlight,
                          alpha=highlight,
                          size=highlight),
  ) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values=c('darkgray','black'),
                     breaks = c(0,1)) +
  scale_alpha_manual(values=c(0.07,0.7),
                     breaks = c(0,1)) +
  scale_size_manual(values=c(0.2, 0.5),
                    breaks = c(0,1)) +
  coord_cartesian(xlim=c(0, 330), ylim=c(-0.1, 1))
# g1_right

pdf(file=paste0("../Results/Plots/Temp/Sensitivity/",
                "Correlation_curves_R",
                ".pdf"),
    width = 6, height = 2)
print(g1_right)
dev.off()


### Visualize neural activities ####
mat_tmp = kronecker(diag(rep(1,328)), rep(1,240))[1:ncol(reduced.dFF), ]
mat_tmp = mat_tmp/colSums(mat_tmp)
reduced.dFF_ave = reduced.dFF %*% mat_tmp

data_2 = data.frame()
id_tmp = 1:nrow(locs_mat)
reduced.dFF_ave_left = reduced.dFF_ave[id_tmp, ]
for (i in 1:(length(id_tmp))) {
    dFF_tmp = reduced.dFF_ave_left[i, ]
    tmp = bind_cols(dFF=dFF_tmp,
                    time=1:length(dFF_tmp),
                    node_1=i,
                    LR = ifelse(locs_mat[i,2]<0,yes='left',no='right'))
    data_2 = bind_rows(data_2, tmp)
}

g2 = data_2 %>%
  ggplot() +
  geom_line(mapping = aes(x=time, y=dFF, group=node_1, color=LR),
            alpha=ifelse(data_2$node_1 %in% (12:16),
                         yes = 0.8,
                         no = 0.8),
            size=ifelse(data_2$node_1 %in% (12:16),
                        yes = 0.2,
                        no = 0.2)) +
  theme_bw() +
  theme(legend.position = 'none') +
  coord_cartesian(xlim=c(0, 330), ylim=c(-0, 0.15))
# g2

pdf(file=paste0("../Results/Plots/Temp/Sensitivity/",
                "dFF_ave_curves",
                ".pdf"),
    width = 6, height = 2)
print(g2)
dev.off()




# Fit AR model ------------------------------------------------------------

cor_max_mat_left = apply(cor.full.ave_left, c(2,3),max)
nodepairs = which(cor_max_mat_left<0.2 & cor_max_mat_left>0.15, arr.ind = TRUE)


ind_nodepair = 6
node_1 = nodepairs[ind_nodepair,1]
node_2 = nodepairs[ind_nodepair,2]
cor_tmp = cor.full.ave_left[,node_1,node_2]
plot(cor_tmp,type='l')
arima_fit = auto.arima(cor_tmp, seasonal = FALSE, stepwise=FALSE, approximation = FALSE)
arima_fit$coef



# Type I error with various threshold (Figure S.4) -----

set.seed(811)

N_rep = 5000

# Plots: type I error vs rho, where type I error = (# GP with mean zero & max>rho) / (# GP with mean zero).
GP_mat = replicate(N_rep, simulate(arima_fit))
GP_max = apply(GP_mat,2,max)
rho_vec = seq(0,1,0.1)
type_I_error_vec = rho_vec
for (i in 1:length(rho_vec)) {
  rho = rho_vec[i]
  type_I_error_vec[i] = sum(GP_max>rho) / length(GP_max)
}

g3 = ggplot() +
  geom_line(mapping = aes(x=rho_vec, y=type_I_error_vec), 
            alpha=0.7) +
  # geom_point(mapping = aes(x=rho_vec, y=type_I_error_vec),
  #            alpha=0.7) +
  theme_bw() +
  xlab("") +
  ylab("Rate of type I error") +
  scale_x_continuous( breaks = seq(0,1,0.2)) +
  coord_cartesian()
g3

pdf(file=paste0("../Results/Plots/Temp/Sensitivity/",
                "Type_I_error",
                ".pdf"),
    width = 3, height = 2)
print(g3)
dev.off()


# Typical correlation curves and empirical power (Figure S.5) ------

set.seed(13)

cor_mat = apply(cor.full.ave_left, 1, function(A)c(A[upper.tri(A)]))

### Choose the number of representative correlation trajectories
tot_withinss_vec = c()
for (N_clus in 1:10) {
  tmp = kmeans(cor_mat,N_clus,nstart = 3)$tot.withinss
  tot_withinss_vec = c(tot_withinss_vec, tmp)
}
ggplot()+
  geom_line(mapping=aes(x=1:10, y=tot_withinss_vec),alpha=0.7) +
  geom_point(mapping=aes(x=1:10, y=tot_withinss_vec),alpha=0.7) +
  theme_bw()

### Visualization of typical correlation curves
N_clus = 5
centers = kmeans(cor_mat,N_clus,nstart = 3)$centers
order = order(apply(centers,1,max))
centers = centers[order,]
rownames(centers) = c(1:N_clus)
centers = as_tibble(t(centers))
centers = bind_cols(centers, time=1:nrow(centers))

data_4 = pivot_longer(centers, cols = c(1:N_clus), 
                      names_to = 'Curve', values_to = 'Correlation') %>%
  mutate(Curve = as.factor(as.numeric(Curve)))
g4 = data_4 %>%
  ggplot() +
  geom_line(mapping = aes(x=time, y=Correlation, group=Curve, 
                          color=Curve), alpha=0.8) +
  scale_color_manual(breaks = c(1:N_clus),
                     values = c('darkgray',hue_pal()(N_clus-1))) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlab("Time (min)") +
  ylab("Correlation") +
  scale_y_continuous( breaks = seq(0,1,0.2)) +
  coord_cartesian()
g4

pdf(file=paste0("../Results/Plots/Temp/Sensitivity/",
                "Typical_correlation",
                ".pdf"),
    width = 3, height = 2.5)
print(g4)
dev.off()





# Obtain and visualize empirical power under various threshold
rho_vec = seq(0.1-0.1,0.9+0.1,0.1)
data_5 = c()
for (k in 1:N_clus) {
  center = unlist(centers[,k])
  GP_mat = replicate(N_rep, center+simulate(arima_fit))
  GP_max = apply(GP_mat,2,max)
  power_vec = rho_vec
  for (i in 1:length(rho_vec)) {
    rho = rho_vec[i]
    power_vec[i] = sum(GP_max>rho) / length(GP_max)
  }
  tmp = bind_cols(power=power_vec, rho=rho_vec, curve=k)
  data_5 = bind_rows(data_5, tmp)
}

g5 = data_5 %>%
  filter(curve > 1) %>%
  mutate(curve = as.factor(curve)) %>%
  ggplot() +
  geom_line(mapping = aes(x=rho, y=power, group=curve, color=curve),
            alpha=0.7) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlab("") +
  ylab("Power") +
  coord_cartesian(xlim = c(0,1))

pdf(file=paste0("../Results/Plots/Temp/Sensitivity/",
                "Power",
                ".pdf"),
    width = 3, height = 2.5)
print(g5)
dev.off()
