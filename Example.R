### The default working directory is the current source file location.


# 0. Import all functions and load libraries ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)
library(ggplot2)
library(grDevices)
library(scales)
library(combinat)
library(cluster)
library(mclust)
library(reshape2)


# 1. Generate networks-------------------------------------------------------

### Parameters used for network generation
SEED = 183  # Random seed
N_node = 30 # Number of nodes per subject 
N_clus = 3  # Number of clusters 
total_time = 200  # Total observation time length 
conn_patt_sep = 1.9  # Separability of intensities across clusters (beta)
conn_prob_mean = 0.9  # Connecting probability 
time_shift_mean_vec = rep(40, N_clus) # Mean value of time shifts (corresponding to W/2 in the main text)
t_vec = seq(0,200,length.out=200) # Time grid of connecting intensities 


### Network generation 
network_list = generate_network(SEED = SEED, 
                                    N_node = N_node, 
                                    N_clus = N_clus, 
                                    total_time = total_time, 
                                    conn_patt_sep = conn_patt_sep, 
                                    conn_prob_mean = conn_prob_mean, 
                                    time_shift_mean_vec = time_shift_mean_vec,
                                    t_vec = t_vec)

edge_time_mat_list = network_list$edge_time_mat_list # Generated connecting time matrix

image(edge_time_mat_list[[1]]) # Heatmap of the first generated connecting time matrix 


# 2. Estimation ------------------------------------
gamma = 0.01 # Relative importance of scales of intensities
freq_trun = Inf # Cut-off frequency (corresponding to \ell_0 in the main text). For cumulative-intensity-based algorithm, set freq_trun to be Inf.
N_clus_tmp = 3 # Desired number of clusters 
MaxIter = 10 # Maximal number of iterations between the centering and aligning steps
conv_thres = 0.01 # Threshold in convergence criterion
max_iter = 10 # Maximal number of iterations between updating time shifts and intensities in the centering step


### Apply proposed initialization scheme
res = get_init(edge_time_mat_list = edge_time_mat_list,
                  N_clus = N_clus_tmp,
                  t_vec = t_vec)
clusters_list_init = res$clusters_list # Initial clustering
n0_vec_list_init = res$n0_vec_list # Initial node-specific time shifts
n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init) # Initial edge-specific time shifts


### Apply SidSBM-C 
res = do_cluster_cdf(edge_time_mat_list = edge_time_mat_list, 
                     N_clus = N_clus_tmp, gamma=gamma,
                     total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                     clusters_list_init = clusters_list_init,
                     n0_vec_list_init = n0_vec_list_init, n0_mat_list_init = n0_mat_list_init,
                     freq_trun = freq_trun,
                     conv_thres = conv_thres,
                     MaxIter = MaxIter)

clusters_list_est = res$clusters_list # Estimated clusters
v_vec_list_est = res$v_vec_list # Estimated time shifts
center_cdf_array_est = res$center_cdf_array # Estimated cumulative connecting intensities
center_pdf_array_est = res$center_pdf_array # Estimated connecting intensities


# 3. Visualization ---------------------------------------------------------

### Match the estimated clusters and the true clusters 
### Identifiable upon to permutation
### ONLY for simulation 
cdf_true_array = network_list$cdf_true_array # True cumulative connecting intensities
permn = find_permn(center_cdf_array_from = center_cdf_array_est, 
                   center_cdf_array_to = cdf_true_array)$permn


### Plot estimated intensities (dashed curves) and true intensities (solid curves)
### In tabular legend, "L" represent truth, "R" represent estimation.
center_pdf_array_est_permn = center_pdf_array_est[permn, permn, ] # Estimated intensities with permuted cluster label
pdf_true_array = network_list$pdf_true_array # True connecting intensities
g = plot_pdf_array(pdf_array_list = list(pdf_true_array,
                                         center_pdf_array_est_permn), 
                   clus_size_vec_1 = rep(N_node/N_clus, N_clus),
                   clus_size_vec_2 = sapply(clusters_list_est[[1]][permn], length),
                   t_vec = t_vec, 
                   subj1_name = 'Truth', subj2_name = 'Est',
                   y_lim = c(0,max(c(center_pdf_array_est_permn, pdf_true_array))) )$g
grid.arrange(g)

### Plot estimated time shifts vs true time shifts
v_true_list = network_list$time_shift_list # True time shifts
plot(unlist(v_true_list), unlist(v_vec_list_est), 
     xlab = "True time shifts", ylab = "Estimated time shifts"); 
abline(a=0, b=1, col=2)


# 4. Evaluate estimation error -----------------------------------------------

### Compute MISE of connecting intensities (i.e., $f$)
dist_mat = matrix(nrow=N_clus, ncol=N_clus)
for (q in 1:N_clus) {
  for (k in 1:N_clus) {
    dist_mat[q,k] = sqrt(sum( (center_pdf_array_est_permn[q,k,] - pdf_true_array[q,k,])^2 * (t_vec[2]-t_vec[1]) ))
  }
}
f_mean_sq_err = mean(dist_mat[upper.tri(dist_mat, diag=TRUE)]^2)
print(paste("Mean Integrated Squared Error (MISE) of estimated intensities:", f_mean_sq_err))

### Compute ARI: accuracy of clustering (i.e., $z$) 
membership_true_list = network_list$membership_true_list # True cluster memberships
ARI_tmp = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est[[1]]), 
                      memb_true_vec = membership_true_list[[1]])
print(paste("Adjusted Rand Index (ARI) of estimated memberships:", ARI_tmp))


# 5. Perform model selection ---------------------------------------------------------
gamma = 0.001
N_clus_min = N_clus - 1 # Minimum candidate number of clusters (i.e., $K$)
N_clus_max = N_clus + 1 # Maximum candidate number of clusters (i.e., $K$)
freq_trun_vec = c(2,3,4) # Candidate values of frequency truncation (i.e., $\ell_0$).
res_list = list() # List of fitted results for candidates of $K$ and $\ell_0$ 

### Fit models for candidate numbers of clusters and candidate frequency truncations
for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
  res_list[[ind_N_clus]] = list()
  N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
  for (ind_freq_trun in 1:length(freq_trun_vec)) {
    freq_trun = freq_trun_vec[ind_freq_trun]
    ### Apply proposed initialization scheme
    res = get_init(edge_time_mat_list = edge_time_mat_list,
                      N_clus = N_clus_tmp,
                      t_vec = t_vec)
    clusters_list_init = res$clusters_list # Initial clustering
    n0_vec_list_init = res$n0_vec_list # Initial node-specific time shifts
    n0_mat_list_init = n0_vec2mat(n0_vec = n0_vec_list_init) # Initial edge-specific time shifts
    ### Apply SidSBM-P 
    res = do_cluster_pdf(edge_time_mat_list = edge_time_mat_list, 
                         N_clus = N_clus_tmp, gamma = gamma,
                         total_time = total_time, max_iter=max_iter, t_vec=t_vec,
                         clusters_list_init = clusters_list_init,
                         n0_vec_list_init = n0_vec_list_init, n0_mat_list_init = n0_mat_list_init,
                         freq_trun = freq_trun,
                         conv_thres = conv_thres,
                         MaxIter = MaxIter)
    ### Save fitted results for current $K$ and $\ell_0$
    res_list[[ind_N_clus]][[ind_freq_trun]] = res 
  }
}

### Conduct model selection using ICL
sel_mod_res = select_model(edge_time_mat_list = edge_time_mat_list,
                           N_node_vec = N_node,
                           N_clus_min = N_clus_min,
                           N_clus_max = N_clus_max,
                           result_list = res_list,
                           total_time = total_time)
### Print the estimated number of clusters
ind_best_N_clus = sel_mod_res$ind_best_N_clus
N_clus_est = c(N_clus_min:N_clus_max)[ind_best_N_clus] 
print(paste("Estimated number of clusters:",N_clus_est))

### Print the estimated frequency truncation parameter
ind_best_freq_trun = sel_mod_res$ind_best_freq_trun
freq_trun_est = freq_trun_vec[ind_best_freq_trun] 
print(paste("Estimated frequency truncation:", freq_trun_est))

### Print the ICL values of all candidate models
ICL_mat = sel_mod_res$ICL_mat
colnames(ICL_mat) = paste("N_clus =",c(N_clus_min:N_clus_max))
rownames(ICL_mat) = paste("freq_trun =", freq_trun_vec)
print(ICL_mat) 


