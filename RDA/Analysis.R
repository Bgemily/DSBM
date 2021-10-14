
# Default working directory: Sources/

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(tidyverse)
library(scales)
library(combinat)
library(cluster)
library(mclust)
library(reshape2)

# Get network information for the L/R side ------------------------------------


data_folder = "../Processed_FunctionalData/"
path_vec = list.files(data_folder, full.names = TRUE)
file_vec = list.files(data_folder, full.names = FALSE)
path_vec = path_vec[1]

edge_time_mat_list = vector(mode = "list", length = length(path_vec)*2)


for(m in 1:length(path_vec)){ 
  path = path_vec[m]
  
  ### Read information from data
  edge_time_mat_L = as.matrix(read.csv(paste(path, '/EdgeTime_L.csv', sep='')))
  edge_time_mat_L = edge_time_mat_L[,-1]
  edge_time_mat_R = as.matrix(read.csv(paste(path, '/EdgeTime_R.csv', sep='')))
  edge_time_mat_R = edge_time_mat_R[,-1]

  edge_time_mat_list[c(2*m-1,2*m)] = list(edge_time_mat_L, edge_time_mat_R)
}





# Apply algorithm ---------------------------------------------------------


N_clus = 3 # Number of clusters
MaxIter = 10 # Maximal iteration number 
bw = 5 # Smoothing bandwidth
freq_trun = 5 # Cut-off frequency


max_time = max(sapply(edge_time_mat_list, function(edge_time_mat)max(edge_time_mat[which(edge_time_mat<Inf)])))
t_vec = seq(0, max_time+10, 1)

clusters_list = center_pdf_array_list = v_vec_list = vector("list",length(edge_time_mat_list))

for (subj in 1:length(edge_time_mat_list)) {
  res = get_init_v4(edge_time_mat_list = edge_time_mat_list[subj], 
                    N_clus = N_clus, 
                    t_vec = t_vec)
  
  clusters_list_init = res$clusters_list
  n0_vec_list_init = res$n0_vec_list
  v_vec_list_init = res$v_vec_list
  n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init)
  
  res = do_cluster_v14(edge_time_mat_list = edge_time_mat_list[subj], 
                       N_clus = N_clus, 
                       clusters_list_init = clusters_list_init, 
                       n0_vec_list_init = n0_vec_list_init, 
                       n0_mat_list_init = n0_mat_list_init, 
                       freq_trun=freq_trun, MaxIter = MaxIter,
                       t_vec = t_vec, bw = bw)
  
  clusters_list[subj] = res$clusters_list
  v_vec_list[subj] = res$v_vec_list
  
  ### Get estimated pdf using kernel smoothing
  v_mat_list = n0_vec2mat(n0_vec = res$v_vec_list)
  n0_mat_list = lapply(v_mat_list, function(v_mat)round(v_mat/(t_vec[2]-t_vec[1])))
  center_pdf_array = get_center_pdf_array_v2(edge_time_mat_list = edge_time_mat_list[subj], 
                                             clusters_list = res$clusters_list, 
                                             n0_mat_list = n0_mat_list, 
                                             t_vec = t_vec, bw=bw)
  center_pdf_array_list[[subj]] = center_pdf_array
  
  ### Save result
  folder_path = paste0('../Results/Rdata/real_data_results/', file_vec[(subj+1)%/%2])
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  file_name = ifelse(subj%%2==1, yes = "Left", no = "Right") 
  now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
  save(list(clusters_list = clusters_list[[subj]],
            center_pdf_array_list = center_pdf_array_list[[subj]],
            v_vec_list = v_vec_list[[subj]]), 
       file = paste0(folder_path, '/', file_name, '_', now_trial, '.Rdata'))
  
}



# Visualization -----------------------------------------------------------


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


