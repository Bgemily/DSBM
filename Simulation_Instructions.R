# 1. Instructions for reproducing results in Figure 5 and Figure 6 -----------------------------------
### To reproduce the simulation results for SidSBM-C, run a for loop within which one of the following variables
### N_node (i.e., $p$), conn_patt_sep (i.e., $beta$), and time_shift_mean (i.e., $W/2$)
### is allowed to change and the other two are fixed. Other parameters are the same as the above example.

### The variables, when allowed to change, take values in the following sets:
# N_node = {30, 42, 54, 66, 78, 90} 
# conn_patt_sep = {1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9} 
# time_shift_mean_vec = {rep(40,N_clus), rep(35,N_clus), rep(30,N_clus), rep(25,N_clus), rep(20,N_clus), rep(15,N_clus)}

### Default values of variables:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1.3
# time_shift_mean_vec = rep(40,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# t_vec = seq(0,200,length.out=200) 

### To reproduce the simulation results for SidSBM-P, use function 'do_cluster_pdf()' instead of 'do_cluster_cdf()'.



# 2. Instructions for reproducing the results in Figure 7(a) -----------------------------------
### For random initialization, use get_init_random() instead of get_init().
### When calling function do_cluster_cdf(), set argument save_est_history = TRUE.



# 3. Instructions for reproducing the results in Figure 7(b) -----------------------------------
### Apply do_cluster_cdf() to simulated data with N_clus = 1, 2, 3, 4, 5, and freq_trun = Inf
### Apply do_cluster_pdf() to simulated data with N_clus = 1, 2, 3, 4, 5, and freq_trun = 2, 3, 4
### Apply select_model() with N_clus_min = 1, N_clus_max = 5


# Instructions for reproducing the results in Figure S1 and S2 -------------------------
### For Scenario A, set the arguments of generate_network() as follows:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1.9
# time_shift_mean_vec = rep(40,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# t_vec = seq(0,200,length.out=200) 

### For Scenario B, set the arguments of generate_network() as follows:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1.9
# time_shift_mean_vec = rep(40,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# conn_prob_mean = 0.5
# conn_prob_rad = 0.4
# t_vec = seq(0,200,length.out=200) 

### For Scenario C, set the arguments of generate_network() as follows:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1
# time_shift_mean_vec = rep(40,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# conn_prob_mean = 0.5
# conn_prob_rad = 0.4
# t_vec = seq(0,200,length.out=200) 

### Apply do_cluster_cdf() with the following varying arguments:
# gamma = 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10
# N_clus = 2, 3, 4


# Instructions for reproducing the results in Figure S3 -------------------------
### Set the arguments of generate_network() as follows:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1.9
# time_shift_mean_vec = rep(40,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# t_vec = seq(0,200,length.out=200) 

### To reproduce the results in Figure S3(a), apply do_cluster_cdf() with the following arguments:
# gamma = 0.001, 0.01, 0.1, 1, 10
# N_clus = 3

### To reproduce the results in Figure S3(b), apply do_cluster_pdf() with the following arguments:
# gamma = 10^(-6), 10^(-5), 10^(-4), 10^(-3), 10^(-2)
# N_clus = 3
# freq_trun = 4


# Instructions for reproducing the results in Figure S4 -------------------------
### Set the arguments of generate_network() as follows:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1.9
# time_shift_mean_vec = rep(40,N_clus), rep(35,N_clus), rep(30,N_clus), rep(25,N_clus), rep(20,N_clus), rep(15,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# t_vec = seq(0,200,length.out=200) 

### Apply do_cluster_cdf() with the following arguments:
# gamma = 0.01
# N_clus = 3

### Apply do_cluster_pdf() with the following arguments:
# gamma = 0.0001
# N_clus = 3
# freq_trun = 4


# Instructions for reproducing the results in Figure S5 -------------------------
### Set the arguments of generate_network() as follows:
# N_node = 30, 90
# N_clus = 3
# conn_patt_sep = 1.9
# time_shift_mean_vec = rep(40,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# t_vec = seq(0,200,length.out=200) 

### Apply do_cluster_cdf() with the following arguments:
# gamma = 0.01, 0.03, 0.1, 0.3, 0.4, ..., 1.3, 0.93, 0.97,  
# N_clus = 3


# Instructions for reproducing the results in Figure S6 -------------------------
### To reproduce results of the shape invariant model (SIM), 
### apply get_init() and do_cluster_cdf() with N_clus = 1



# Instructions for reproducing the results in Figure S7 and S8 -------------------------
### To reproduce the simulation results, run a for loop within which one of the 
### following variables is allowed to change and the other one is fixed.
###  N_node (i.e., $p$) and conn_patt_sep (i.e., $beta$)

### The two variables, when allowed to change, take values in the following sets:
# N_node = {30, 42, 54, 66, 78, 90} 
# conn_patt_sep = {1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9} 

### Default values of variables are as follows:
# N_node = 30
# N_clus = 3
# conn_patt_sep = 1.3
# time_shift_mean_vec = rep(0,N_clus)
# total_time = 200  
# conn_prob_mean = 0.9  
# t_vec = seq(0,200,length.out=200) 

### In the step of initialization, use get_init_jittertruev() instead of get_init().
### Call do_cluster_cdf() to apply SidSBM-C, and do_clusters_pdf() to apply SidSBM-P.
### When calling do_cluster_cdf() or do_cluster_pdf(), set argument fix_timeshift = TRUE.



