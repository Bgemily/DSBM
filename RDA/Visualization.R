#!/usr/bin/env Rscript



# Import functions --------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Visualization (NeuIPS submission) -----------------------------------------------------------

### Read in data analysis results
data_folder = "../Results/Rdata/RDA/do_cluster_v14.2.1/"
path_vec = list.files(data_folder, full.names = TRUE, recursive = TRUE)

res_list = vector(mode = "list", length = length(path_vec))

for(m in 1:length(path_vec)){ 
  path = path_vec[m]
  ### Read in results from data
  load(path)
  res_list[[m]] = res
}


### Collect from all subjects: clusters_list, center_pdf_array
clusters_list = lapply(res_list, function(res) res$clusters_list)
center_pdf_array_list = lapply(res_list, function(res) res$center_pdf_array)


### Align clusters and intensities
clusters_list[[1]] = clusters_list[[1]][c(2,3,1)]
clusters_list[[2]] = clusters_list[[2]][c(2,1,3)]
center_pdf_array_list[[1]] = center_pdf_array_list[[1]][c(2,3,1),c(2,3,1),]
center_pdf_array_list[[2]] = center_pdf_array_list[[2]][c(2,1,3),c(2,1,3),]

N_clus = length(clusters_list[[1]])
for (q in 1:N_clus) {
  for (k in 1:N_clus) {
    center_pdf_array_list[[1]][q,k,] = shift_v2(center_pdf_array_list[[1]][q,k,], n0 = -17)
  }
}


### Visualize estimated connecting intensities (Figure 3 in DynamicNetworks.pdf) ----
clus_size_vec = rowSums(sapply(clusters_list, function(clusters)sapply(clusters,length)))
print(plot_pdf_array_v2(pdf_array_list = center_pdf_array_list, 
                        pdf_true_array = (center_pdf_array_list[[1]]+center_pdf_array_list[[2]])/2,
                        clus_size_vec = clus_size_vec,
                        t_vec = t_vec, y_lim = c(0,max(unlist(center_pdf_array_list)))))




# Visualization -----------------------------------------------------------

max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
t_vec = seq(0, max_time+10, length.out = 1000)


dir.create(paste0('../Results/Plots/Temp/Real_data_summary/'), recursive = TRUE, showWarnings = FALSE)

### Align clusters and intensities
clusters_list[[1]] = clusters_list[[1]][c(2,3,1)]
clusters_list[[2]] = clusters_list[[2]][c(2,1,3)]
center_pdf_array_list[[1]] = center_pdf_array_list[[1]][c(2,3,1),c(2,3,1),]
center_pdf_array_list[[2]] = center_pdf_array_list[[2]][c(2,1,3),c(2,1,3),]

N_clus = length(clusters_list[[1]])
for (q in 1:N_clus) {
  for (k in 1:N_clus) {
    center_pdf_array_list[[1]][q,k,] = shift_v2(center_pdf_array_list[[1]][q,k,], n0 = -17)
  }
}


### F: Connection pattern  ----
pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'F','.pdf'), width = 3.5, height = 3.5)
clus_size_vec = rowSums(sapply(clusters_list, function(clusters)sapply(clusters,length)))
print(plot_pdf_array_v2(pdf_array_list = center_pdf_array_list, 
                        pdf_true_array = (center_pdf_array_list[[1]]+center_pdf_array_list[[2]])/2,
                        clus_size_vec = clus_size_vec,
                        t_vec = t_vec, y_lim = c(0,max(unlist(center_pdf_array_list)))))
dev.off()

### v: time shift  ----
pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'v_est_vs_init','.pdf'), width = 3, height = 3)
plot(unlist(v_vec_list_init[avai_subj][1]),unlist(v_vec_list[1]))
abline(a=0,b=1,col=2)
dev.off()

pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'v_distribution','.pdf'), width = 2, height = 1.5)
membership_list = lapply(clusters_list[1], clus2mem)

data = tibble(v_vec=unlist(v_vec_list[1]), 
              membership=unlist(membership_list), 
              subj=rep(1:length(v_vec_list[1]), times=sapply(v_vec_list[1],length)))
# g<- data %>%
#   ggplot(aes(x=v_vec,
#                      group=(membership),
#                      color=as.factor(membership)))+
#   # geom_histogram( aes(), color="#e9ecef", alpha=1, position = 'dodge', bins=5 ) +
#   # geom_boxplot(aes(color=as.factor(membership)), alpha=1) +
#   geom_density( aes(), alpha=0.6, position = 'dodge', bw=20) +
#   scale_fill_manual(values=palette()[1+1:3]) +
#   theme_bw()+
#   theme(legend.position = "none") +
#   # facet_wrap(~subj) +
#   # scale_y_continuous(position="right") +
#   theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
#   xlab('Activation time (v)')
# print(g)

g<- data %>%
  ggplot(aes(x=v_vec,group=as.factor(membership),
                     color=as.factor(membership)))+
  geom_density( aes(), alpha=0.6, position = 'dodge', bw=20) +
  scale_fill_manual(values=palette()[1+1:3]) +
  theme_bw()+
  theme(legend.position = "none") +
  # facet_wrap(~subj) +
  scale_x_continuous(limits = c(0,250)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  xlab('Activation time (v)')

print(g)
dev.off()


### Connection time matrix  ----
pdf(paste0("../Results/Plots/Temp/Real_data_summary/",'Edge_time_noise_no_order','.pdf'),width = 4,height = 4)
v_mat_list = n0_vec2mat(v_vec_list)
plot_edge_time_mat(edge_time_mat = edge_time_mat_list[[1]], clusters = clusters_list[[1]],
                   n0_mat = v_mat_list[[1]]/(t_vec[2]-t_vec[1]), t_vec = t_vec, 
                   zlim = c(0,max(unlist(edge_time_mat_list)[unlist(edge_time_mat_list)<Inf])),
                   denoise = FALSE, reorder = FALSE)
dev.off()


pdf(paste0("../Results/Plots/Temp/","Real_data_summary/",'Edge_time_denoise_reorder','.pdf'),width = 4,height = 4)
v_mat_list = n0_vec2mat(v_vec_list)
plot_edge_time_mat(edge_time_mat = edge_time_mat_list[[1]], clusters = clusters_list[[1]],
                   n0_mat = v_mat_list[[1]]/(t_vec[2]-t_vec[1]), t_vec = t_vec, 
                   zlim = c(0,max(unlist(edge_time_mat_list)[unlist(edge_time_mat_list)<Inf])),
                   denoise = TRUE, reorder = TRUE)
dev.off()

pdf(paste0("../Results/Plots/Temp/","Real_data_summary/",'Edge_time_denoise_no_order','.pdf'),width = 4,height = 4)
v_mat_list = n0_vec2mat(v_vec_list)
plot_edge_time_mat(edge_time_mat = edge_time_mat_list[[1]], clusters = clusters_list[[1]],
                   n0_mat = v_mat_list[[1]]/(t_vec[2]-t_vec[1]), t_vec = t_vec, 
                   zlim = c(0,max(unlist(edge_time_mat_list)[unlist(edge_time_mat_list)<Inf])),
                   denoise = TRUE, reorder = FALSE)
dev.off()

pdf(paste0("../Results/Plots/Temp/","Real_data_summary/",'Edge_time_noise_reorder','.pdf'),width = 4,height = 4)
v_mat_list = n0_vec2mat(v_vec_list)
plot_edge_time_mat(edge_time_mat = edge_time_mat_list[[1]], clusters = clusters_list[[1]],
                   n0_mat = v_mat_list[[1]]/(t_vec[2]-t_vec[1]), t_vec = t_vec, 
                   zlim = c(0,max(unlist(edge_time_mat_list)[unlist(edge_time_mat_list)<Inf])),
                   denoise = FALSE, reorder = TRUE)
dev.off()



pdf(paste0("../Results/Plots/Temp/Real_data_summary/",'Color_bar','.pdf'),width = 3,height = 4)
colorBreaks = c(seq(1,20,length.out=10),seq(21,100,length.out=300),seq(101,300,length.out=300))
fields::image.plot(edge_time_mat_list[[1]], zlim=c(0,max(unlist(edge_time_mat_list)[unlist(edge_time_mat_list)<Inf])), 
                   breaks = colorBreaks,
                   col=fields::tim.colors(length(colorBreaks)-1),
                   horizontal = TRUE,
                   xaxt="n",yaxt="n")
dev.off()


### Spatial location  ----

for (i in 1:length(membership_list)) {
  m = avai_subj[i]
  N_node = nrow(locs_mat_list[[m]])
  mem_tmp = rep(0,N_node)
  mem_tmp[non_iso_inds_list[[m]]] = membership_list[[i]]
  membership_list[[i]] = mem_tmp
}

data = as.data.frame(reduce(locs_mat_list[avai_subj][1:1],rbind))
colnames(data) = c('AP','LR','DV')
data = cbind(data, 
             membership=as.factor(unlist(membership_list)),
             subj=rep(1:length(membership_list), times=sapply(membership_list,length)))
data = pivot_longer(data, cols = c('AP','LR','DV'), names_to = "coord", values_to = "loc")

pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'spatial_location','.pdf'), 
    width = 2, height = 4)
# lims = data.frame(coord=c("AP","DV","LR"), xmin=c(-0.5, 0, -2),xmax=c(9,140,3))
g<- ggplot(data, aes(x=loc,
                     group=membership,
                     color=membership))+
  geom_density( aes(), alpha=0.6, position = 'dodge') +
  scale_color_manual(values=palette()[c(8,2:4)]) +
  theme_bw()+
  theme(legend.position = "none") +
  facet_wrap(~coord, scales='free',ncol=1) +
  xlab('') +
  # scale_x_continuous(limits = c(-2, 6))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  scale_y_continuous(position = "right")
print(g)
dev.off()



### mnx decomposition  ----

data = tibble(subj=rep(1:length(membership_list[1]),sapply(membership_list[1],length)),
              membership=unlist(membership_list[1]),
              mnx=as.factor(unlist(mnx_vec_list[avai_subj][1])), islet=NA) 
for (m in 1:length(avai_subj[1])) {
  if (!is.null(islet_vec_list[avai_subj][[m]])) {
    data[which(data$subj==m), 'islet'] = as.numeric(islet_vec_list[avai_subj][[m]])
  }
}

pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'mnx_decomp','.pdf'), 
    width = 2, height = 4)
data_2 = data %>% 
  mutate(mnx=factor(mnx,levels=c(1,0)), islet=as.factor(as.numeric(islet)), 
         membership=factor(membership, levels=c(0,3,2,1))) 
data_2 %>%
  ggplot(mapping = aes(fill=membership,x=mnx)) +
  # ggplot(mapping = aes(fill=membership,x=factor(1))) +
  scale_fill_manual(values=palette()[c(8,4,3,2)]) +
  geom_bar(position="fill") +
  # facet_wrap(facets=. ~ mnx, ncol = 1,
  #            labeller = labeller(mnx = c(`0`=paste0('mnx- (p=', as.numeric(data_2 %>% filter(mnx==0) %>% summarise(n())), ')'), 
  #                                        `1`=paste0('mnx+ (p=', as.numeric(data_2 %>% filter(mnx==1) %>% summarise(n())), ')')))) +
  # coord_polar(theta="y") +
  geom_text(data = data_2 %>% 
              count(mnx, membership) %>% group_by(mnx) %>% 
              summarise(perc=round(n/sum(n),3), membership=membership) %>% 
              arrange(desc(membership), .by_group = TRUE) ,
            aes(x = mnx, label = scales::percent(perc), y=perc, fill=NULL), 
            position = position_fill(0.5), size=2.5)+
  scale_x_discrete(labels=c("mnx+","mnx-"))+
  xlab("")+
  theme_bw()+
  theme(legend.position = 'none', 
        panel.grid  = element_blank(),
        # axis.text = element_blank(),
        axis.text.y = element_blank(), 
        axis.title = element_blank(),  
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank())

dev.off()


### cluster composition using mnx+/-  ####
pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'clus_decomp','.pdf'), 
    width = 4, height = 2)
data_2 = data %>%
  mutate(mnx=as.factor(mnx), islet=as.factor(as.numeric(islet)), membership=as.factor(membership))
data_2 %>%
  # ggplot(mapping = aes(fill=membership,x=mnx)) +
  ggplot(mapping = aes(fill=mnx,x=factor(1))) +
  scale_fill_manual(values=palette()[c(8,2:6)]) +
  geom_bar(position="fill") +
  facet_grid(facets=. ~ membership,
             labeller = labeller(membership = c(`0`=paste0('Iso (p=', as.numeric(data_2 %>% filter(membership==0) %>% summarise(n())), ')'), 
                                         `1`=paste0('Red (p=', as.numeric(data_2 %>% filter(membership==1) %>% summarise(n())), ')'),
                                         `2`=paste0('Green (p=', as.numeric(data_2 %>% filter(membership==2) %>% summarise(n())), ')'),
                                         `3`=paste0('Blue (p=', as.numeric(data_2 %>% filter(membership==3) %>% summarise(n())), ')')))) +
  coord_polar(theta="y") +
  geom_text(data = data_2 %>% 
              count(membership, mnx) %>% group_by(membership) %>% 
              summarise(perc=round(n/sum(n),2), mnx=mnx) %>% 
              arrange(desc(mnx), .by_group = TRUE) ,
            aes(x = 1.1, label = scales::percent(perc), y=perc, fill=NULL), 
            position = position_fill(0.5), size=2.5)+
  xlab("")+
  theme_bw()+
  theme(legend.position = 'none', 
        panel.grid  = element_blank(),
        axis.text = element_blank(),
        axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks.y = element_blank())

dev.off()



### islet decomposition ----

pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'islet_decomp','.pdf'), 
    width = 3, height = 3)
data_3 = data %>% 
  mutate(mnx=as.factor(mnx), islet=as.factor((islet*1)), membership=as.factor(membership)) %>%
  filter(mnx!=0 | islet==0 | is.na(islet)) %>% 
  filter(subj %in% which(sapply(islet_vec_list[avai_subj][1:2],length)>0) ) 
data_3 %>%
  ggplot(mapping = aes(fill=membership,x=islet)) +
  scale_fill_manual(values=palette()[c(8,2:6)]) +
  geom_bar(position="fill") +
  geom_text(data = data_3 %>% count(islet, membership) %>% group_by(islet) %>% 
              summarise(perc=round(n/sum(n),2), membership=membership) %>% 
              arrange(desc(membership), .by_group = TRUE) ,
            aes(x = islet, label = scales::percent(perc), y=perc, fill=NULL), 
            position = position_fill(0.5), size=2.5)+
  geom_text(data = data_3 %>% group_by(islet) %>% summarise(n()), 
            aes(x = islet, label =  `n()`, y = 1, fill = NULL), vjust=-0.2, size=2.5)+
  scale_x_discrete(labels=c("Ventral IN","MN","CoLA/CoSA"))+
  theme_bw()+
  theme(legend.position = 'none', 
        # axis.title.x = element_text(size=2), axis.text.x = element_text(size=2), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
dev.off()



### dFF curves ----
avai.inds = as.matrix(data.table::fread(paste('../Processed_data/func_20150410/AvaiNeurons.csv',sep='')))
avai.inds=avai.inds[,-1];
dat.dFF = as.matrix(data.table::fread(paste('../Processed_data/func_20150410/dFF.csv',sep='')))
dat.dFF=dat.dFF[,-1]
reduced.dFF=dat.dFF[avai.inds,]
locs.all = as.matrix(read.csv(paste('../Processed_data/func_20150410/locs_all.csv',sep='')))
locs.all = locs.all[,-1]
locs=locs.all[avai.inds,]

g_tmp = plot_local_traces(locs=locs_mat_list[[1]], membership=membership_list[[1]],
                           reduced.dFF=reduced.dFF[which(locs[,2]<0),], 
                           snapshot_mins=c(30,120,240), window_length=1,scale=4)
g_tmp
g_tmp = gridExtra::arrangeGrob(grobs = g_tmp, ncol = length(g_tmp))
ggsave(filename = paste0('../Results/Plots/Temp/Real_data_summary/Traces.pdf'), plot = g_tmp,
       width = 6, height = 3.5, units = "in")



### Network snapshots ----
tmp = rep(1, nrow(edge_time_mat_list[[1]]))
node_pair=c(22,23)
tmp[node_pair] = RColorBrewer::brewer.pal(3,name="Dark2")[1:2]
cols = t(col2rgb(tmp))
vertex.size = rep(3, nrow(edge_time_mat_list[[1]]))
vertex.size[node_pair] = 6
v_mat_list = n0_vec2mat(n0_vec = v_vec_list)

plot_network_animation(locs = locs_mat_list[[1]][non_iso_inds_list[[1]],], 
             edge.time = edge_time_mat_list[[1]]-v_mat_list[[1]],
             output = paste0('../Results/Plots/Temp/Real_data_summary/'), 
             # filename = paste0(avai_subj), 
             vertex.size = vertex.size,
             window_list = list(c(0,48),c(0,58),c(0,68)), 
             asp=1, save_plots = T, delay=120,
             alpha = 150,
             cols = cols)

### Network animation ----
tmp = rep(1, nrow(edge_time_mat_list[[1]]))
node_pair=c(22,23)
tmp[node_pair] = RColorBrewer::brewer.pal(3,name="Dark2")[1:2]
cols = t(col2rgb(tmp))
vertex.size = rep(3, nrow(edge_time_mat_list[[1]]))
# vertex.size[node_pair] = 6

v_mat_list = n0_vec2mat(n0_vec = v_vec_list)

plot_network_animation(locs = locs_mat_list[[1]][non_iso_inds_list[[1]],], 
                       edge.time = edge_time_mat_list[[1]]-v_mat_list[[1]],
                       output = paste0('../Results/Plots/Temp/Real_data_summary/'), 
                       # filename = paste0(avai_subj), 
                       vertex.size = vertex.size,
                       window_list = lapply(c(46:50),function(x)c(0,x)), 
                       asp=1, save_plots = T, delay=80,
                       alpha = 800,
                       cols = cols*0,
                       )


### plot edge times as lollipop plots ----

# node_pair=c(24,26)
node_pair=c(22,23)
# node_pair=sample(clusters_list[[1]][[1]],2)

### V2
  

library(tidyverse)
edge_time_tibble = as_tibble(t(edge_time_mat_list[[1]][node_pair, ]))
edge_time_tibble = pivot_longer(edge_time_tibble, 
                                        cols=everything(), 
                                        names_to = 'node_id', 
                                        values_to = 'edge_time') %>%
                      mutate(node_id=factor(node_id, levels=c('V1','V2')))

colorBreaks = c(seq(1,20,length.out=10),seq(21,100,length.out=300),seq(101,300,length.out=300))
edge_time_rescale = sapply((edge_time_tibble %>% 
                              filter(edge_time<Inf))$edge_time, 
                           function(edge_time) sum(colorBreaks<edge_time)/length(colorBreaks)*300)

g_raw = edge_time_tibble %>% 
  filter(edge_time<Inf) %>%
  mutate(y=1, edge_time=jitter(edge_time, amout=5)) %>%
  ggplot() +
  # geom_segment( aes(x=edge_time, xend=edge_time, y=0, yend=y, 
  #                   group=node_id, color=node_id, linetype=factor(node_id,levels=c('V2','V1'))),
  #               size=0.5) +
  geom_segment( aes(x=edge_time, xend=edge_time, y=0, yend=y),
                color=fields::tim.colors(300)[edge_time_rescale],
                size=0.5) +
  facet_wrap(~node_id, nrow=2)+
  xlab('Connecting time (min)') +
  # theme_void()+
  xlim(0,300) +
  scale_color_brewer(palette="Dark2", direction = -1) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'none')
ggsave('../Results/Plots/Temp/Real_data_summary/edge_time_spike_train.pdf',
       plot = g_raw,
       width = 3,height = 1.5)



shifted_edge_time_tibble = as_tibble(t(round(edge_time_mat_list[[1]][node_pair, ] - 
                                               v_mat_list[[1]][node_pair, ])))
shifted_edge_time_tibble = pivot_longer(shifted_edge_time_tibble, 
                                        cols=everything(), 
                                        names_to = 'node_id', 
                                        values_to = 'edge_time') %>%
                              mutate(node_id=as.factor(node_id))
colorBreaks = c(seq(1,20,length.out=10),seq(21,100,length.out=300),seq(101,300,length.out=300))
edge_time_rescale = sapply((shifted_edge_time_tibble %>% 
                              filter(edge_time<Inf))$edge_time, 
                           function(edge_time) sum(colorBreaks<edge_time)/length(colorBreaks)*300)

g_shift = shifted_edge_time_tibble %>% 
  filter(edge_time<Inf) %>%
  mutate(y=1, edge_time=jitter(edge_time,amount = 5)) %>%
  ggplot() +
  # geom_segment( aes(x=edge_time, xend=edge_time, y=0, yend=y, 
  #                   group=node_id, color=node_id, linetype=node_id),
  #               size=0.5) +
  geom_segment( aes(x=edge_time, xend=edge_time, y=0, yend=y),
                color=fields::tim.colors(300)[edge_time_rescale],
                size=0.5) +
  facet_wrap(~node_id, nrow=2)+
  xlab('Time post maturation (min)') +
  xlim(c(0,300)) +
  scale_color_brewer(palette="Dark2",) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'none')
  
ggsave('../Results/Plots/Temp/Real_data_summary/edge_time_spike_train_shifted.pdf',
       plot = g_shift,
       width = 3,height = 1.5)




### plot cluster membership matrix -----

pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'mem_raw_order','.pdf'), 
    width = 6, height = 2)
image(as.matrix(membership_list[[1]][1:44]), col=gray.colors(3,start=1,end=0), xaxt="n",yaxt="n")
dev.off()

pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'mem_reorder','.pdf'), 
    width = 6, height = 2)
image(as.matrix(sort(membership_list[[1]][1:44])), col=gray.colors(3,start=1,end=0), xaxt="n",yaxt="n")
dev.off()

mem_mat = dummies::dummy(membership_list[[1]][1:44])
image(t(mem_mat), col=gray.colors(2,start = 1,end=0))


# Choose N_clus for real data ---------------------------------------------------------------


ICL_vec = c()
file_list = list.files(path = '../Results/Rdata/real_data_results/',
                       pattern = 'Partial_subj_N_clus*', full.names = TRUE)
for (file in file_list) {
  load(file)
  ICL = get_ICL(edge_time_mat_list = edge_time_mat_list[avai_subj],
                clusters_list = clusters_list,
                v_vec_list = v_vec_list,
                center_pdf_array = center_pdf_array)
  ICL_vec = c(ICL_vec, ICL)
}
plot(2:5,ICL_vec, type='b')




############################################



