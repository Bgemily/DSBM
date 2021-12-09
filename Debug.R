# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# Read in data analysis results ----------------------

### Read in data analysis results
data_folder = paste0(params$data_folder, params$subj_name, "/") 
path_vec = list.files(data_folder, full.names = TRUE, recursive = TRUE)

res_list = vector(mode = "list", length = length(path_vec))

for(m in 1:length(path_vec)){ 
  path = path_vec[m]
  ### Read in results from data
  load(path)
  res_list[[m]] = res
}


### Collect from all subjects: clusters_list, center_pdf_array, v_vec
clusters_list = lapply(res_list, function(res) res$clusters_list)
center_pdf_array_list = lapply(res_list, function(res) res$center_pdf_array)
v_vec_list = lapply(res_list, function(res) res$v_vec)
t_vec = res$t_vec
n0_mat_list = lapply(v_vec_list, function(v_vec) n0_vec2mat(n0_vec = v_vec/(t_vec[2]-t_vec[1])))



### Get connecting probabilities ----
order_clus_1 = order(rowSums(apply(center_pdf_array_list[[1]], 
                                   MARGIN = c(1,2),
                                   function(vec)sum(vec*(t_vec[2]-t_vec[1])))), 
                     decreasing = TRUE)
order_clus_2 = order(rowSums(apply(center_pdf_array_list[[2]], 
                                   MARGIN = c(1,2),
                                   function(vec)sum(vec*(t_vec[2]-t_vec[1])))), 
                     decreasing = TRUE)

### Visualize estimated connecting intensities ------
clus_size_list = lapply(clusters_list, function(clusters)sapply(clusters,length))

pdf_array_1 = center_pdf_array_list[[1]][order_clus_1,order_clus_1, ,drop=FALSE]#[1:3,1:3,]
pdf_array_2 = center_pdf_array_list[[2]][order_clus_2,order_clus_2, ,drop=FALSE]#[1:3,1:3,]



g1 = plot_pdf_array_v2(pdf_array_list = center_pdf_array_list[[1]][order_clus_1,order_clus_1, ,drop=FALSE], 
                       pdf_true_array = center_pdf_array_list[[1]][order_clus_1,order_clus_1, ,drop=FALSE],
                       clus_size_vec = clus_size_list[[1]][order_clus_1],
                       t_vec = t_vec, x_lim = c(), y_lim = c(0,max(unlist(center_pdf_array_list[[1]]))))

g2 = plot_pdf_array_v2(pdf_array_list = center_pdf_array_list[[2]][order_clus_2,order_clus_2, ,drop=FALSE], 
                       pdf_true_array = center_pdf_array_list[[2]][order_clus_2,order_clus_2, ,drop=FALSE],
                       clus_size_vec = clus_size_list[[2]][order_clus_2],
                       t_vec = t_vec, x_lim = c(), y_lim = c(0,max(unlist(center_pdf_array_list[[2]]))))

grid.arrange(g1, g2, ncol=2)
