#!/usr/bin/env Rscript


rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Estimate and save edge_time_mat ------------------------------------------------------------------

# data_folder = "../Processed_data/"
# path_vec = list.files(data_folder, full.names = TRUE)
# 
# window_length = 240
# rho = 0.4
# 
# for (path in path_vec) {
#   rm(cor.full.ave)
#   load(paste0(path,'/cor_full_ave','_win',window_length,'.rdata'))
#   
#   adj.full = cor.full.ave;
#   for(i in 1:(dim(cor.full.ave)[1])){
#     adj.full[i,,] = cor.full.ave[i,,]>rho;
#   }
#   
#   
#   ### Find the earliest connecting time for all pairs of neurons 
#   edge.time = matrix(Inf, nrow=dim(cor.full.ave)[2], ncol=dim(cor.full.ave)[3])
#   for(i in 1:nrow(edge.time)){
#     for(j in 1:ncol(edge.time)){
#       tmp = min(which(adj.full[,i,j]==1)) # index of interval
#       edge.time[i,j] = ifelse(tmp<Inf, tmp*(window_length/240),Inf) # edge time (using min as the unit)
#     }
#   }
#   diag(edge.time) = Inf
#   
#   write.csv(edge.time, paste(path,'/EdgeTime.csv',sep=''), col.names=F)
#   
#   
# }
# 
# 


# Get network information for the L/R side of all subjects ------------------------------------

data_folder = "../Processed_data/"
path_vec = list.files(data_folder, full.names = TRUE)

file_name_vec = list.files(data_folder)
subj_name_list = sapply(file_name_vec, function(file_name)c(paste0(file_name,'_','L'), 
                                                            paste0(file_name,'_','R')))
subj_name_list = c(subj_name_list)[1]


edge_time_mat_list = locs_mat_list = 
  mnx_vec_list = 
  iso_inds_list = non_iso_inds_list = 
  vector(mode = "list", length = length(subj_name_list))

for(m in 1:1){
  path = path_vec[m]
  
  ### Read information from data
  edge_time_mat = as.matrix(read.csv(paste(path, '/EdgeTime.csv', sep='')))
  edge_time_mat = edge_time_mat[,-1]
  
  avai.inds = as.matrix(read.csv(paste(path,'/AvaiNeurons.csv',sep='')))
  avai.inds = avai.inds[,-1];
  
  locs.all = as.matrix(read.csv(paste(path,'/locs_all.csv',sep='')))
  locs.all = locs.all[,-1]
  locs_mat = locs.all[avai.inds,]
  
  mnx.all = as.matrix(read.csv(paste(path,'/mnx.csv',sep='')))
  mnx.all = mnx.all[,-1]
  mnx_vec = mnx.all[avai.inds]
  
  ### Split each zebrafish into left and right
  inds_L = which(locs_mat[,2]<0)
  inds_R = which(locs_mat[,2]>0)
  
  edge_time_mat_L = edge_time_mat[inds_L, inds_L]
  edge_time_mat_R = edge_time_mat[inds_R, inds_R]
  
  locs_mat_L = locs_mat[inds_L,]
  locs_mat_R = locs_mat[inds_R,]
  
  mnx_vec_L = mnx_vec[inds_L]
  mnx_vec_R = mnx_vec[inds_R]
  
  
  iso_inds_L = which(rowSums(edge_time_mat_L<Inf)==0)
  iso_inds_R = which(rowSums(edge_time_mat_R<Inf)==0)
  
  non_iso_inds_L = setdiff(1:nrow(edge_time_mat_L), iso_inds_L)
  non_iso_inds_R = setdiff(1:nrow(edge_time_mat_R), iso_inds_R)
  
  edge_time_mat_L = edge_time_mat_L[non_iso_inds_L, non_iso_inds_L]
  edge_time_mat_R = edge_time_mat_R[non_iso_inds_R, non_iso_inds_R]
  
  edge_time_mat_list[c(2*m-1,2*m)] = list(edge_time_mat_L, edge_time_mat_R)
  locs_mat_list[c(2*m-1,2*m)] = list(locs_mat_L, locs_mat_R)
  mnx_vec_list[c(2*m-1,2*m)] = list(mnx_vec_L, mnx_vec_R)
  
  iso_inds_list[c(2*m-1,2*m)] = list(iso_inds_L, iso_inds_R)
  non_iso_inds_list[c(2*m-1,2*m)] = list(non_iso_inds_L, non_iso_inds_R)
  
}

### Output of this section: 
### locs_mat_list, mnx_vec_list, islet_vec_list, \\
### edge_time_mat_list, iso_inds_list, non_iso_inds_list

# plot(NULL,xlim=c(0,350),ylim=c(0,0.01))
# for(m in 1:length(edge_time_mat_list)){
#   lines(density(edge_time_mat_list[[m]], bw = 15),col=m)
# }


# Apply algorithm ---------------------------------------------------------


# N_clus = 3
# MaxIter = 10
bw = 5
freq_trun = 5
N_clus = 3

max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
t_vec = seq(0, max_time+10, length.out = 1000)
t_vec = seq(0, max_time+10, 1)

clusters_list = center_pdf_array_list = v_mat_list = vector("list",2)

for (avai_subj in 1:2) {
  res = get_init_v4(edge_time_mat_list = edge_time_mat_list[avai_subj], 
                    N_clus = N_clus, 
                    t_vec = t_vec)
  
  clusters_list_init = res$clusters_list
  n0_vec_list_init = res$n0_vec_list
  v_vec_list_init = res$v_vec_list
  n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
  
  res = do_cluster_v14(edge_time_mat_list = edge_time_mat_list[avai_subj], 
                       N_clus = N_clus, 
                       clusters_list_init = clusters_list_init, 
                       n0_vec_list_init = n0_vec_list_init, 
                       n0_mat_list_init = n0_mat_list_init, 
                       freq_trun=freq_trun, MaxIter = 10,
                       t_vec = t_vec, bw = bw)
  
  clusters_list[avai_subj] = res$clusters_list
  
  # 
  # res$clusters_history -> clusters_history
  # res$v_vec_list -> v_vec_list
  # res$center_pdf_array -> center_pdf_array
  # res$center_pdf_array_list -> center_pdf_array_list
  
  
  ### Get estimated pdf using kernel smoothing
  v_mat_list_tmp = n0_vec2mat(n0_vec = res$v_vec_list)
  n0_mat_list_tmp = lapply(v_mat_list_tmp, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
  center_pdf_array = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list[avai_subj], 
                                             clusters_list = res$clusters_list, 
                                             n0_mat_list = n0_mat_list_tmp, 
                                             t_vec = t_vec, bw=bw)
  center_pdf_array_list[[avai_subj]] = center_pdf_array
  v_mat_list[avai_subj] = v_mat_list_tmp
  
  
}

### Save results
# save(list=c('subj_name_list', 'avai_subj',
#             "edge_time_mat_list", 'iso_inds_list', 'non_iso_inds_list',
#             'locs_mat_list','mnx_vec_list',
#             'clusters_list', 'v_vec_list_init',
#             'center_pdf_array', 'center_pdf_array_list'),
#      file=paste0('../Results/Rdata/real_data_results/',
#                  'Single_subj_two_networks', '_', 'N_clus', N_clus, '.Rdata'))







# Visualization -----------------------------------------------------------

max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
t_vec = seq(0, max_time+10, length.out = 1000)


dir.create(paste0('../Results/Plots/Temp/Real_data_summary/'), recursive = TRUE, showWarnings = FALSE)

### NeurIPS submission ----
pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'F','.pdf'), width = 3.5, height = 3.5)
clusters_list[[1]] = clusters_list[[1]][c(2,3,1)]
clusters_list[[2]] = clusters_list[[2]][c(2,1,3)]
center_pdf_array_list[[1]] = center_pdf_array_list[[1]][c(2,3,1),c(2,3,1),]
center_pdf_array_list[[2]] = center_pdf_array_list[[2]][c(2,1,3),c(2,1,3),]

### Manually align 
for (q in 1:N_clus) {
  for (k in 1:N_clus) {
    center_pdf_array_list[[1]][q,k,] = shift(center_pdf_array_list[[1]][q,k,], n0 = 17)
  }
}

clus_size_vec = rowSums(sapply(clusters_list, function(clusters)sapply(clusters,length)))
pdf(file = paste0("../Results/Plots/Temp/Real_data_summary/", 'F','.pdf'), width = 3.5, height = 3.5)
print(plot_pdf_array_v2(pdf_array_list = center_pdf_array_list, 
                        pdf_true_array = (center_pdf_array_list[[1]]+center_pdf_array_list[[2]])/2,
                        clus_size_vec = clus_size_vec,
                        t_vec = t_vec, y_lim = c(0,max(unlist(center_pdf_array_list)))))
dev.off()

membership_list = lapply(clusters_list, clus2mem)
for (m in 1:2) {
  N_node = nrow(locs_mat_list[[m]])
  mem_tmp = rep(0,N_node)
  mem_tmp[non_iso_inds_list[[m]]] = membership_list[[i]]
  membership_list[[i]] = mem_tmp
}
contingency_table = table(unlist(membership_list), unlist(mnx_vec_list))
colSums(contingency_table)


