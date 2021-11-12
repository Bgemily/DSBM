### Based on v1
### Compatible with new simulation result structure (i.e., nested folders)
### Remove redundant visualization
### 

# Import functions --------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library("tidyverse")
library("dplyr")
library(ggplot2)


# cdf + jitter (V==0) ------------------------------------------------------------
path_vec = rep(0,5)

path_vec[1] = "../Results/Rdata/SNR_Vis0/main_v5_cdf_timeshift_jitter_v4/jitter_time_rad/0/pr=1,n=30,beta=1.2/"
path_vec[2] = "../Results/Rdata/SNR_Vis0/main_v5_cdf_timeshift_jitter_v4/jitter_time_rad/5/pr=1,n=30,beta=1.2/"
path_vec[3] = "../Results/Rdata/SNR_Vis0/main_v5_cdf_timeshift_jitter_v4/jitter_time_rad/10/pr=1,n=30,beta=1.2/"
path_vec[4] = "../Results/Rdata/SNR_Vis0/main_v5_cdf_timeshift_jitter_v4/jitter_time_rad/15/pr=1,n=30,beta=1.2/"
path_vec[5] = "../Results/Rdata/SNR_Vis0/main_v5_cdf_timeshift_jitter_v4/jitter_time_rad/20/pr=1,n=30,beta=1.2/"

param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err", "ICL_vec")))
  results_df = bind_rows(bind_cols(results_list[[1]],"jitter_level"="0"),
                         bind_cols(results_list[[2]],"jitter_level"="5"),
                         bind_cols(results_list[[3]],"jitter_level"="10"),
                         bind_cols(results_list[[4]],"jitter_level"="15"),
                         bind_cols(results_list[[5]],"jitter_level"="20"))
  
  ### Manipulate column "param_value"
  results_df = results_df %>% 
    mutate(param_value = factor(param_value, levels = sort(unique(param_value), 
                                                           decreasing = param_name=="V")),
           jitter_level = as_factor(jitter_level)) 
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse","V_mse")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'vs', '_', "N_node", ".pdf"), 
        width = 4, height = 4)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "V_mse" = v_mean_sq_err), 
                 color=jitter_level)) +
      stat_summary(aes(group=jitter_level), position = position_dodge(.3),
                   geom="pointrange",
                   fun.min = function(x)quantile(x,0.25),
                   fun.max = function(x)quantile(x,0.75),
                   fun = median) +
      stat_summary(aes(group=jitter_level),position = position_dodge(.3),
                   geom="line",
                   fun = "median") +
      theme(legend.position = "bottom") +
      # guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      scale_y_continuous(limits = switch(measurement,
                                         "1-ARI" = c(0,1),
                                         "f_mse" = c())) +
      ylab(measurement) +
      xlab("N_node")
    
    print(g)
    dev.off()
  }
  
}



# cdf vs pdf (V!=0) ------------------------------------------------------------
path_vec = rep(0,5)

path_vec[1] = "../Results/Rdata/SNR_Vnot0/main_v5_cdf/pr=1,n=30,beta=1.2/"
path_vec[2] = "../Results/Rdata/SNR_Vnot0/main_v5_pdf_v3_select_freqtrun/pr=1,n=30,beta=1.2/"
path_vec[3] = "../Results/Rdata/SNR_Vnot0/main_v5_pdf_v4_pairwise_alignment/pr=1,n=30,beta=1.2/"

param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err", "ICL_vec")))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="CDF"), 
                         bind_cols(results_list[[2]],"method"="PDF+simult_align"), 
                         bind_cols(results_list[[3]],"method"="PDF+pairwise_align"))
  
  ### Manipulate column "param_value"
  results_df = results_df %>% 
    mutate(param_value = factor(param_value, levels = sort(unique(param_value), 
                                                           decreasing = param_name=="V")) ) 
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse","V_mse")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'vs', '_', "N_node", ".pdf"), 
        width = 4, height = 4)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "V_mse" = v_mean_sq_err), 
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.3),
                   geom="pointrange",
                   fun.min = function(x)quantile(x,0.25),
                   fun.max = function(x)quantile(x,0.75),
                   fun = median) +
      stat_summary(aes(group=method),position = position_dodge(.3),
                   geom="line",
                   fun = "median") +
      theme(legend.position = "bottom") +
      # guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      scale_y_continuous(limits = switch(measurement,
                                         "1-ARI" = c(0,1),
                                         "f_mse" = c())) +
      ylab(measurement) +
      xlab("N_node")
    
    print(g)
    dev.off()
  }
  
}



# pdf performance vs N_node (V!=0) ------------------------------------------------------------
path_vec = rep(0,5)

path_vec[1] = "../Results/Rdata/SNR_Vnot0/main_v5_pdf_v4_pairwise_alignment/pr=1,n=30,beta=1.2/"

param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err", "ICL_vec")))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="PDF+pairwise_align"))
  
  ### Manipulate column "param_value"
  results_df = results_df %>% 
    mutate(param_value = factor(param_value, levels = sort(unique(param_value), 
                                                           decreasing = param_name=="V")) ) 
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse","V_mse")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'vs', '_', "N_node", ".pdf"), 
        width = 4, height = 4)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "V_mse" = v_mean_sq_err), 
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.3),
                   geom="pointrange",
                   fun.min = function(x)quantile(x,0.25),
                   fun.max = function(x)quantile(x,0.75),
                   fun = median) +
      stat_summary(aes(group=method),position = position_dodge(.3),
                   geom="line",
                   fun = "median") +
      theme(legend.position = "bottom") +
      # guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      scale_y_continuous(limits = switch(measurement,
                                         "1-ARI" = c(0,1),
                                         "f_mse" = c())) +
      ylab(measurement) +
      xlab("N_node")
    
    print(g)
    dev.off()
  }
  
}



# ICL/log_lik/penalty vs N_clus and freq_trun -----------------------------

path_vec = rep(0,1)

path_vec[1] = "../Results/Rdata/SNR_Vnot0/main_v5_v4_multifreqtrun/pr=0.4,n=30,beta=1.3/"
param_name_vec = list.files(path_vec[1])
param_name_vec = c("n/90/freqtrun")

val_n = 90 # 30, 42, 54, 66, 78, 90
val_beta = 1.8
### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: freq_trun | ICL | log-lik | penalty
  func_tmp = function(folder_path, param_name) {
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name),
                           measurement = c("ICL_vec","compl_log_lik_vec", "penalty_vec", "penalty_2_vec"))
  }
  results_list = lapply(path_vec, func_tmp, param_name=param_name)
  results_df = results_list[[1]] %>% 
    pivot_longer(cols = starts_with("ICL_vec"), names_to = "N_clus_ICL", values_to = "ICL") %>% 
    pivot_longer(cols = starts_with("compl_log_lik_vec"), names_to = "N_clus_log_lik", values_to = "log_lik") %>% 
    pivot_longer(cols = starts_with("penalty_vec"), names_to = "N_clus_penalty", values_to = "penalty") %>%
    pivot_longer(cols = starts_with("penalty_2_vec"), names_to = "N_clus_penalty_2", values_to = "penalty_2")
  
  ### Plot ICL vs N_clus
  for (measurement in c("ICL","log_lik","penalty","penalty_2")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    switch(param_name, "beta/1.8/freqtrun"=paste0("Beta_",1.8), 
                           "n/90/freqtrun"=paste0("N_node_",90), 
                           "V"="V",""), '_', 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement), ".pdf"), 
        width = 4, height = 2.5)
    g = results_df %>% 
      mutate(N_basis = 2*param_value+1) %>%
      mutate(N_basis = as.factor(N_basis)) %>%
      ggplot(aes(x=switch(measurement,
                          "ICL" = N_clus_ICL,
                          "log_lik" = N_clus_log_lik,
                          "penalty_2" = N_clus_penalty_2,
                          "penalty" = N_clus_penalty), 
                 y=switch(measurement,
                          "ICL" = ICL,
                          "log_lik" = log_lik,
                          "penalty_2" = penalty_2,
                          "penalty" = penalty), 
                 color=N_basis)) +
      stat_summary(aes(group=N_basis), position = position_dodge(.2),
                   geom="pointrange",
                   fun = mean,
                   fun.min = function(x) mean(x)-sd(x),
                   fun.max = function(x) mean(x)+sd(x) ) +
      stat_summary(aes(group=N_basis),position = position_dodge(.2),
                   geom="line",
                   fun = "mean") +
      # theme(legend.position = "none") +
      scale_y_continuous() +
      ylab(measurement) +
      scale_x_discrete(labels=1:5) +
      xlab("Number of clusters")
    
    print(g)
    dev.off()
  }
}


# ICL/log_like/penalty vs N_clus -----------------------------------------

path_vec = rep(0,2)

path_vec[1] = "../Results/Rdata/SNR_Vis0/apply_ppsbm_ICL/pr=0.4,n=30,beta=1.3/"
path_vec[2] = "../Results/Rdata/SNR_Vis0/main_v5_v5_adap_freq/pr=0.4,n=30,beta=1.3/"


param_name_vec = list.files(path_vec[1])
# param_name_vec = c("beta","n")

val_n = 90 # 30, 42, 54, 66, 78, 90
val_beta = 1.8

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | correct_N_clus | method
  func_tmp = function(folder_path, param_name) {
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name),
                           measurement = c("ICL_vec","compl_log_lik_vec", "penalty_vec"))
  }
  results_list = lapply(path_vec, func_tmp, param_name=param_name)
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="ppsbm"),
                         bind_cols(results_list[[2]],"method"="main_v5"),) %>% 
    filter(param_value == switch(param_name, "beta"=val_beta, "n"=val_n)) %>%
    pivot_longer(cols = starts_with("ICL_vec"), names_to = "N_clus_ICL", values_to = "ICL") %>% 
    pivot_longer(cols = starts_with("compl_log_lik_vec"), names_to = "N_clus_log_lik", values_to = "log_lik") %>% 
    pivot_longer(cols = starts_with("penalty_vec"), names_to = "N_clus_penalty", values_to = "penalty")
  
  ### Plot ICL vs N_clus
  for (measurement in c("ICL","log_lik","penalty")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    switch(param_name, "beta"=paste0("Beta_",val_beta), 
                           "n"=paste0("N_node_",val_n), 
                           "V"="V"), '_', 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement), ".pdf"), 
        width = 4, height = 2.5)
    g = results_df %>% 
      ggplot(aes(x=switch(measurement,
                          "ICL" = N_clus_ICL,
                          "log_lik" = N_clus_log_lik,
                          "penalty" = N_clus_penalty), 
                 y=switch(measurement,
                          "ICL" = ICL,
                          "log_lik" = log_lik,
                          "penalty" = penalty), 
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.2),
                   geom="pointrange",
                   fun = mean,
                   fun.min = function(x) mean(x)-sd(x),
                   fun.max = function(x) mean(x)+sd(x) ) +
      stat_summary(aes(group=method),position = position_dodge(.2),
                   geom="line",
                   fun = "mean") +
      # theme(legend.position = "none") +
      scale_y_continuous() +
      ylab(measurement) +
      scale_x_discrete(labels=1:5) +
      xlab("Number of clusters")
    
    print(g)
    dev.off()
  }
}


# Cluster number selection (V==0) -----------------------------------------

path_vec = rep(0,1)

path_vec[1] = "../Results/Rdata/SNR_Vis0/apply_ppsbm_ICL/pr=0.4,n=30,beta=1.3/"
path_vec[2] = "../Results/Rdata/SNR_Vis0/main_v5_v7_largefreqtrun/pr=0.4,n=30,beta=1.3/"


param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | N_clus_est | method
  results_list = lapply(path_vec, function(folder_path) extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name),
                                                                               measurement = "N_clus_est"))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="ppsbm"),
                         bind_cols(results_list[[2]],"method"="our"),)
  
  ### Manipulate column "param_value"
  results_df = results_df %>% 
    mutate(param_value = factor(param_value, levels = sort(unique(param_value), 
                                                           decreasing = param_name=="V")) ) 
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("N_clus_est")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    switch(param_name, "beta"="Beta", "n"="N_node", "V"="V"), '_', 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement), ".pdf"), 
        width = 4, height = 2.5)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "incorrect_N_clus" = 1-correct_N_clus,
                          "N_clus_est" = N_clus_est,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "V_mse" = v_mean_sq_err), 
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.2),
                   geom="pointrange",
                   fun = mean,
                   fun.min = function(x) mean(x)-sd(x),
                   fun.max = function(x) mean(x)+sd(x) ) +
      stat_summary(aes(group=method),position = position_dodge(.2),
                   geom="line",
                   fun = "mean") +
      # theme(legend.position = "none") +
      scale_y_continuous(limits = switch(measurement,
                                         # "incorrect_N_clus" = c(0,1),
                                         "1-ARI" = c(0,1),
                                         "f_mse" = c(0.0,0.006))) +
      ylab(measurement) +
      xlab(ifelse(param_name=="n", yes="p", no=param_name))
    
    print(g)
    dev.off()
  }
}


results_df %>% 
  filter(param_value==90, method=='ppsbm') %>% 
  summarise(mean=mean(N_clus_est), sd=sd(N_clus_est))
results_df %>% 
  filter(param_value==90, method=='our') %>% 
  summarise(mean=mean(N_clus_est), sd=sd(N_clus_est))

# cdf vs pdf (V!=0) ------------------------------------------------------------


# path_vec = rep(0,2)
# 
# path_vec[1] = "../Results/Rdata/SNR_Vnot0/our_v3.4/p=0.7,n=30,beta=1.5,V=80"
# path_vec[2] = "../Results/Rdata/SNR_Vnot0/our_v3.1.1.1/p=0.7,n=30,beta=1.5,V=80"

path_vec = rep(0,6)

path_vec[1] = "../Results/Rdata/SNR_Vnot0/our_v3.4/p=0.7,n=30,beta=1.5,V=80"
path_vec[2] = "../Results/Rdata/SNR_Vnot0/our_v3.1.1.1/p=0.7,n=30,beta=1.5,V=80"
path_vec[3] = "../Results/Rdata/SNR_Vnot0/our_v3.4.5_rad0/p=0.7,n=30,beta=1.5,V=80/"
path_vec[4] = "../Results/Rdata/SNR_Vnot0/our_v3.4.5_rad25/p=0.7,n=30,beta=1.5,V=80/"
path_vec[5] = "../Results/Rdata/SNR_Vnot0/our_v3.4.5_rad50/p=0.7,n=30,beta=1.5,V=80/"
path_vec[6] = "../Results/Rdata/SNR_Vnot0/our_v3.4.5_rad75/p=0.7,n=30,beta=1.5,V=80/"

param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name)))
  
  # results_df = bind_rows(bind_cols(results_list[[1]],"method"="pdf"), 
  #                             bind_cols(results_list[[2]],"method"="cdf"))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="pdf"), 
                         bind_cols(results_list[[2]],"method"="cdf"),
                         bind_cols(results_list[[3]],"method"="cdf+pdf_rad0"),
                         bind_cols(results_list[[4]],"method"="cdf+pdf_rad25"),
                         bind_cols(results_list[[5]],"method"="cdf+pdf_rad50"),
                         bind_cols(results_list[[6]],"method"="cdf+pdf_rad75"),)
  
  ### Manipulate column "param_value"
  results_df = results_df %>% 
    mutate(param_value = factor(param_value, levels = sort(unique(param_value), 
                                                           decreasing = param_name=="V")) ) 
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse","V_mse")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    switch(param_name, "beta"="Beta", "n"="N_node", "V"="V"), '_', 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement), ".pdf"), 
        width = 6, height = 3)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "V_mse" = v_mean_sq_err), 
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.2),
                   geom="point",
                   fun = mean,
                   fun.min = function(x)quantile(x,0.25),
                   fun.max = function(x)quantile(x,0.75)) +
      stat_summary(aes(group=method),position = position_dodge(.2),
                   geom="line",
                   fun = "mean") +
      coord_cartesian(ylim = switch(measurement,
                                     "1-ARI" = c(0,0.4),
                                     "f_mse" = c(0.0,0.02))) +
      ylab(measurement) +
      xlab(ifelse(param_name=="n", yes="p", no=param_name))
    
    print(g)
    dev.off()
  }
  
  
  
  
}




# cdf vs pdf (V==0) ------------------------------------------------------------


path_vec = rep(0,2)

path_vec[1] = "../Results/Rdata/SNR_Vis0/main_v5_v2/pr=0.4,n=30,beta=1.3"
path_vec[2] = "../Results/Rdata/SNR_Vis0/main_v5_v3_adap_freq/pr=0.4,n=30,beta=1.3"


param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name)))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="Fixed_freq_trun"), 
                         bind_cols(results_list[[2]],"method"="Adaptive_freq_trun"))
  
  ### Manipulate column "param_value"
  results_df = results_df %>% 
    mutate(param_value = factor(param_value, levels = sort(unique(param_value), 
                                                           decreasing = param_name=="V")) ) 
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    switch(param_name, "beta"="Beta", "n"="N_node", "V"="V"), '_', 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement), ".pdf"), 
        width = 4, height = 2.5)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "V_mse" = v_mean_sq_err), 
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.2),
                   geom="pointrange",
                   fun = mean,
                   fun.min = function(x)quantile(x,0.25),
                   fun.max = function(x)quantile(x,0.75)) +
      stat_summary(aes(group=method),position = position_dodge(.2),
                   geom="line",
                   fun = "mean") +
      # theme(legend.position = "none") +
      scale_y_continuous(limits = switch(measurement,
                                         "1-ARI" = c(0,1),
                                         "f_mse" = c(0.0,0.006))) +
      ylab(measurement) +
      xlab(ifelse(param_name=="n", yes="p", no=param_name))
    
    print(g)
    dev.off()
  }
  
  
  
  
}






# Multiple subject ------------------------------------------------------------

library("tidyverse")

path_vec = rep(0,6)

path_vec[1] = "../Results/Rdata/Multi_subj_case1_v3.1.1/n=30,beta=1.5,p=0.8//" # p=0.8
path_vec[2] = "../Results/Rdata/Multi_subj_case2_v3.1.1/n=30,beta=1.5,p=0.8//"
path_vec[3] = "../Results/Rdata/Multi_subj_case3_v3.1.1/n=30,beta=1.5,p=0.8//"
path_vec[4] = "../Results/Rdata/Multi_subj_case4_v3.3.1/n=30,beta=1.5,p=0.8//"


# path_vec[4] = "../Results/Rdata/20210318_120018"
# path_vec[5] = "../Results/Rdata/20210318_120039"
# path_vec[6] = "../Results/Rdata/20210318_120053"

# path_vec[1] = "../Results/Rdata/20210317_151254"
# path_vec[2] = "../Results/Rdata/20210317_151423"
# path_vec[3] = "../Results/Rdata/20210317_151456"


path_vec[c(1,2,3)] = path_vec[c(1,3,2)]

file_list_mat = lapply(path_vec, function(path)list.files(path=path, full.names = T, recursive = TRUE))

  Lambda_mean_sq_err_longlong = data.frame()
  ARI_mean_longlong = data.frame()
  F_mean_sq_err_longlong = data.frame()
  v_mean_sq_err_longlong = data.frame()
  for (j in 1:4) {
    file_list = file_list_mat[[j]]
    
    ARI_mean = extract_measurement(file_list = file_list, measurement = "ARI_mean")[[1]]
    ARI_mean = rename_with(ARI_mean, function(s)as.numeric(gsub("[^0-9.]", "", s)))
    ARI_mean_long = ARI_mean %>% 
      pivot_longer(cols = everything(), 
                   names_to = "N_subj", 
                   values_to = "ARI") %>%
      mutate(Setup = j)
    ARI_mean_longlong = rbind(ARI_mean_longlong, ARI_mean_long)
    
    F_mean_sq_err = extract_measurement(file_list = file_list, measurement = "F_shape_mean_sq_err")[[1]]
    F_mean_sq_err = rename_with(F_mean_sq_err, function(s)as.numeric(gsub("[^0-9.]", "", s)))
    F_mean_sq_err_long = F_mean_sq_err %>% 
      pivot_longer(cols = everything(), 
                   names_to = "N_subj", 
                   values_to = "F_mse") %>%
      mutate(Setup = j)
    F_mean_sq_err_longlong = rbind(F_mean_sq_err_longlong, F_mean_sq_err_long)
    
  }
  
  pdf(file=paste0("../Results/Plots/Temp/", "ARI_vs_N_subj.pdf"), width = 6, height = 4)
  g = ARI_mean_longlong %>%
    # filter(as.numeric(Jitter_radius)<=20) %>%
    mutate(Setup=as.factor(Setup), N_subj=as.factor(as.numeric(N_subj))) %>%
    # mutate(Setup=recode(Setup, 
    #                     `1`="S1,v8", 
    #                     `2`="S2,v8", 
    #                     `3`="S3,v8",
    #                     `4`="S1,v7",
    #                     `5`="S2,v7",
    #                     `6`="S3,v7",)) %>%
    ggplot(mapping = aes(x=N_subj, y=1-ARI, fill=Setup, color=Setup)) +
    # geom_boxplot() +
    # scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
    stat_summary(aes(group=Setup,linetype=Setup,shape=Setup), position = position_dodge(.0),
                 geom="point",
                 # fun.min = function(x)quantile(x,0.25),
                 # fun.max = function(x)quantile(x,0.75),
                 fun = mean) +
    stat_summary(aes(group=Setup,linetype=Setup),position = position_dodge(.0),
                 geom="line",
                 fun = "mean") +
    # ylim(c(0,0.3))+
    # scale_y_continuous(limits=c(0,1))+
    labs(linetype = "Scenario",fill = "Scenario", color="Scenario", shape='Scenario')  +
    scale_color_manual(values=c(`1`=scales::hue_pal()(4)[1],
                                `2`=scales::hue_pal()(4)[1],
                                `3`=scales::hue_pal()(4)[1],
                                `4`=scales::hue_pal()(4)[2],
                                `5`=scales::hue_pal()(4)[2],
                                `6`=scales::hue_pal()(4)[2])) 

  
  print(g)
  dev.off()
  
  
  
  
  pdf(file=paste0("../Results/Plots/Temp/", "F_vs_N_subj.pdf"), width = 6, height = 4)
  g = F_mean_sq_err_longlong %>%
    # filter(as.numeric(Jitter_radius)<=20) %>%
    mutate(Setup=as.factor(Setup), N_subj=as.factor(as.numeric(N_subj))) %>%
    ggplot(mapping = aes(x=N_subj, y=F_mse, fill=Setup, color=Setup)) +
    # geom_boxplot() +
    # scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
    stat_summary(aes(group=Setup,linetype=Setup,shape=Setup), position = position_dodge(.0),
                 geom="point",
                 # fun.min = function(x)quantile(x,0.25),
                 # fun.max = function(x)quantile(x,0.75),
                 fun = median) +
    stat_summary(aes(group=Setup,linetype=Setup),position = position_dodge(.0),
                 geom="line",
                 fun = "median")  +
    labs(linetype = "Scenario",fill = "Scenario", color="Scenario", shape='Scenario')  +
    scale_color_manual(values=c(`1`=scales::hue_pal()(4)[1],
                                `2`=scales::hue_pal()(4)[1],
                                `3`=scales::hue_pal()(4)[1],
                                `4`=scales::hue_pal()(4)[2],
                                `5`=scales::hue_pal()(4)[2],
                                `6`=scales::hue_pal()(4)[2])) 
  print(g)
  dev.off()
  







# Unbalanced clusters -----------------------------------------------------

library("tidyverse")

path_vec = rep(0,3)

path_vec[1] = "../Results/Rdata/Unbalanced_clusters_our_v3.2/p=0.5,n=60,beta=1.5,V=0/" 
path_vec[2] = "../Results/Rdata/Unbalanced_clusters_ppsbm/p=0.5,n=60,beta=1.5,V=0/" 


file_list_mat = lapply(path_vec, function(path)list.files(path=path, full.names = TRUE, recursive = TRUE))

rm(results_vec)
results_vec = c("results1", "results2", 'results3')

ARI_mean_biglist = vector("list", 2)
F_mean_sq_err_biglist = vector("list",2)
for (j in 1:2){
  file_list = file_list_mat[[j]]
  ARI_mean_list = extract_measurement(file_list = file_list, 
                                      measurement = "ARI_mean", 
                                      obj_name_vec = results_vec)
  ARI_mean_biglist[[j]] = ARI_mean_list
  F_mean_sq_err_list = extract_measurement(file_list = file_list, 
                                           measurement = "F_shape_mean_sq_err", 
                                           obj_name = results_vec)
  F_mean_sq_err_biglist[[j]] = F_mean_sq_err_list
}

for (j in 1) {
  Lambda_mean_sq_err_longlong = data.frame()
  ARI_mean_longlong = data.frame()
  F_mean_sq_err_longlong = data.frame()
  v_mean_sq_err_longlong = data.frame()
  for (i in 1:length(results_vec)) {
    obj_name = results_vec[i]

    ARI_mean = ARI_mean_biglist[[j]][[i]]
    ARI_mean = rename_with(ARI_mean, function(s)as.numeric(gsub("[^0-9.]", "", s)))
    ARI_mean_long = ARI_mean %>% 
      pivot_longer(cols = everything(), 
                   names_to = "Jitter_radius", 
                   values_to = "ARI") %>%
      mutate(Setup = i)
    ARI_mean_longlong = rbind(ARI_mean_longlong, ARI_mean_long)
    
    F_mean_sq_err = F_mean_sq_err_biglist[[j]][[i]]
    F_mean_sq_err = rename_with(F_mean_sq_err, function(s)as.numeric(gsub("[^0-9.]", "", s)))
    F_mean_sq_err_long = F_mean_sq_err %>% 
      pivot_longer(cols = everything(), 
                   names_to = "Jitter_radius", 
                   values_to = "F_mse") %>%
      mutate(Setup = i)
    F_mean_sq_err_longlong = rbind(F_mean_sq_err_longlong, F_mean_sq_err_long)
    
  }
  
  ### ARI_mean_longlong: Convert x axis variable to a factor
  ARI_mean_longlong = ARI_mean_longlong %>%
    mutate(Setup=as.factor(Setup), Jitter_radius=as.numeric(Jitter_radius)) %>%
    mutate(Setup=recode(Setup, 
                        `1`="Cluster 1", 
                        `2`="Cluster 2", 
                        `3`="Cluster 3"))
  factor_level = unique(ARI_mean_longlong$Jitter_radius)
  
  ### Sort the factor levels as needed
  factor_level = sort(factor_level, decreasing = TRUE)
  
  ARI_mean_longlong = ARI_mean_longlong %>%
    mutate(Jitter_radius=factor(Jitter_radius,levels = factor_level))
  
  ### F_mean_sq_err_longlong: Convert x axis variable to a factor
  F_mean_sq_err_longlong = F_mean_sq_err_longlong %>%
    mutate(Setup=as.factor(Setup), Jitter_radius=as.numeric(Jitter_radius)) %>%
    mutate(Setup=recode(Setup, 
                        `1`="Cluster 1", 
                        `2`="Cluster 2", 
                        `3`="Cluster 3"))
  F_mean_sq_err_longlong = F_mean_sq_err_longlong %>%
    mutate(Jitter_radius=factor(Jitter_radius,levels = factor_level))
  
  
  
  ### Visualization
  pdf(file=paste0("../Results/Plots/Temp/", 'Unbalance_clusters','_','ARI.pdf'), width = 5, height = 3)
  g = ARI_mean_longlong %>%
    ggplot(mapping = aes(x=Jitter_radius, y=1-ARI, fill=Setup, color=Setup)) +
    # geom_boxplot() +
    # scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
    stat_summary(aes(group=Setup), position = position_dodge(.0),
                 geom="point",
                 fun = mean) +
    stat_summary(aes(group=Setup),position = position_dodge(.0),
                 geom="line",
                 fun = "mean") +
    # theme(legend.position = "none") +
    xlab('Largest cluster size') +
    ylim(c(0,1))
  
  print(g)
  dev.off()
  
  
  
  
  pdf(file=paste0("../Results/Plots/Temp/", 'Unbalance_clusters','_', "F.pdf"), width = 5, height = 3)
  g = F_mean_sq_err_longlong %>%
    ggplot(mapping = aes(x=Jitter_radius, y=F_mse, fill=Setup, color=Setup)) +
    # geom_boxplot() +
    # scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
    stat_summary(aes(group=Setup), position = position_dodge(.0),
                 geom="point",
                 fun = mean) +
    stat_summary(aes(group=Setup),position = position_dodge(.0),
                 geom="line",
                 fun = "mean") +
    # theme(legend.position = "none") +
    xlab('Largest cluster size') +
    scale_y_continuous(trans = "log10",limits=c(0.05,200))
  
  print(g)
  dev.off()
  
}



# Jitter_radius ----------------------------------------------------------
library("tidyverse")

path_vec = rep(0,4)
path_vec[1] = "../Results/Rdata/20210323_204751" ### p=0.7
path_vec[2] = "../Results/Rdata/20210323_204837"
path_vec[1] = "../Results/Rdata/20210317_113420" ### p=0.5
path_vec[2] = "../Results/Rdata/20210317_113457"
# 
# path_vec[1] = "../Results/Rdata/20210317_020309" ### p=1
# path_vec[2] = "../Results/Rdata/20210317_020404"


path_vec[1] = "../Results/Rdata/Jitter_radius/Adaptive_order/" 
path_vec[2] = "../Results/Rdata/Jitter_radius/Fixed_order/"


file_list_mat = lapply(path_vec, function(path)list.files(path=path, full.names = T, recursive = TRUE))

Lambda_mean_sq_err_longlong = data.frame()
ARI_mean_longlong = data.frame()
F_mean_sq_err_longlong = data.frame()
v_mean_sq_err_longlong = data.frame()
for (j in 1:2) {
  file_list = file_list_mat[[j]]
  
  ARI_mean_list = extract_measurement(file_list = file_list, measurement = "ARI_mean", obj_name_vec = c('results1'))
  ARI_mean = ARI_mean_list[[1]]
  ARI_mean = rename_with(ARI_mean, function(s)as.numeric(gsub("[^0-9.]", "", s)))
  ARI_mean_long = ARI_mean %>% 
    pivot_longer(cols = everything(), 
                 names_to = "Jitter_radius", 
                 values_to = "ARI") %>%
    mutate(Setup = j)
  ARI_mean_longlong = rbind(ARI_mean_longlong, ARI_mean_long)
  
  F_mean_sq_err = extract_measurement(file_list = file_list, measurement = "F_mean_sq_err", obj_name_vec = c('results1'))
  F_mean_sq_err = F_mean_sq_err[[1]]
  F_mean_sq_err = rename_with(F_mean_sq_err, function(s)as.numeric(gsub("[^0-9.]", "", s)))
  F_mean_sq_err_long = F_mean_sq_err %>% 
    pivot_longer(cols = everything(), 
                 names_to = "Jitter_radius", 
                 values_to = "F_mse") %>%
    mutate(Setup = j)
  F_mean_sq_err_longlong = rbind(F_mean_sq_err_longlong, F_mean_sq_err_long)
  
  v_mean_sq_err = extract_measurement(file_list = file_list, measurement = "v_mean_sq_err", obj_name_vec = c('results1'))
  v_mean_sq_err = v_mean_sq_err[[1]]
  v_mean_sq_err = rename_with(v_mean_sq_err, function(s)as.numeric(gsub("[^0-9.]", "", s)))
  v_mean_sq_err_long = v_mean_sq_err %>% 
    pivot_longer(cols = everything(), 
                 names_to = "Jitter_radius", 
                 values_to = "v_mse") %>%
    mutate(Setup = j)
  v_mean_sq_err_longlong = rbind(v_mean_sq_err_longlong, v_mean_sq_err_long)
}

### for j=3&4: ppsbm ###

# file_list = file_list_mat[[3]]
# 
# ARI_mean_tmp = extract_measurement(file_list = file_list, measurement = "ARI_mean", obj_name = 'results6')
# ARI_mean_tmp = rename_with(ARI_mean_tmp, function(s)as.numeric(gsub("[^0-9.]", "", s)))
# ARI_mean_tmp = ARI_mean_tmp$`1`
# 
# F_mean_sq_err_tmp = extract_measurement(file_list = file_list, measurement = "F_mean_sq_err", obj_name = 'results6')
# F_mean_sq_err_tmp = rename_with(F_mean_sq_err_tmp, function(s)as.numeric(gsub("[^0-9.]", "", s)))
# F_mean_sq_err_tmp = F_mean_sq_err_tmp$`1`
# 
# file_list = file_list_mat[[4]]
# 
# ARI_mean_tmp_2 = extract_measurement(file_list = file_list, measurement = "ARI_mean", obj_name = 'results6')
# ARI_mean_tmp_2 = rename_with(ARI_mean_tmp_2, function(s)as.numeric(gsub("[^0-9.]", "", s)))
# ARI_mean_tmp_2 = ARI_mean_tmp_2$`1`
# ARI_mean_tmp = c(ARI_mean_tmp, ARI_mean_tmp_2)
# 
# F_mean_sq_err_tmp_2 = extract_measurement(file_list = file_list, measurement = "F_mean_sq_err", obj_name = 'results6')
# F_mean_sq_err_tmp_2 = rename_with(F_mean_sq_err_tmp_2, function(s)as.numeric(gsub("[^0-9.]", "", s)))
# F_mean_sq_err_tmp_2 = F_mean_sq_err_tmp_2$`1`
# F_mean_sq_err_tmp = c(F_mean_sq_err_tmp, F_mean_sq_err_tmp_2)
# 
# ARI_mean_2 = ARI_mean
# for (j in 1:ncol(ARI_mean)) {
#   ARI_mean_2[,j] = sample(ARI_mean_tmp, size=nrow(ARI_mean), replace = TRUE)
# }
# ARI_mean_long = ARI_mean_2 %>% 
#   pivot_longer(cols = everything(), 
#                names_to = "Jitter_radius", 
#                values_to = "ARI") %>%
#   mutate(Setup = 3)
# ARI_mean_longlong = rbind(ARI_mean_longlong, ARI_mean_long)
# 
# 
# F_mean_sq_err_2 = F_mean_sq_err
# for (j in 1:ncol(F_mean_sq_err)) {
#   F_mean_sq_err_2[,j] = sample(F_mean_sq_err_tmp, size=nrow(F_mean_sq_err), replace = TRUE)
# }
# F_mean_sq_err_long = F_mean_sq_err_2 %>% 
#   pivot_longer(cols = everything(), 
#                names_to = "Jitter_radius", 
#                values_to = "F_mse") %>%
#   mutate(Setup = 3)
# F_mean_sq_err_longlong = rbind(F_mean_sq_err_longlong, F_mean_sq_err_long)

# ARI_mean_longlong = filter(ARI_mean_longlong, as.numeric(Jitter_radius)%%4==2)
# F_mean_sq_err_longlong = filter(F_mean_sq_err_longlong, as.numeric(Jitter_radius)%%4==2)
# v_mean_sq_err_longlong = filter(v_mean_sq_err_longlong, as.numeric(Jitter_radius)%%4==2)


pdf(file="../Results/Plots/Temp/ARI_vs_jitter_rad.pdf", width = 6, height = 4)
ARI_mean_longlong %>%
  # filter(as.numeric(Jitter_radius)<=20) %>%
  mutate(Setup=as.factor(Setup), Jitter_radius=as.factor(as.numeric(Jitter_radius))) %>%
  mutate(Jitter_radius=factor(Jitter_radius,levels=rev(levels(Jitter_radius)))) %>%
  mutate(Setup=recode(Setup, 
                      `1`="Adaptive order", 
                      `3`="Adaptive order", 
                      `2`="Fixed order",
                      `4`="Fixed order",)) %>%
  mutate(Setup=factor(Setup,levels=c('Adaptive order','PPSBM','Fixed order'))) %>%
  ggplot(mapping = aes(x=Jitter_radius, y=1-ARI, fill=Setup, color=Setup)) +
  stat_summary(aes(group=Setup), position = position_dodge(.0),
               geom="point",
               fun = mean) +
  stat_summary(aes(group=Setup),position = position_dodge(.0),
               geom="line",
               fun = "mean") +
  scale_color_manual(values=c('Adaptive order'=scales::hue_pal()(4)[1],'Fixed order'=scales::hue_pal()(4)[2])) +
  xlab("Jitter level") 
  # ylim(c(0,1))
  


dev.off()




pdf(file="../Results/Plots/Temp/F_vs_jitter_rad.pdf", width = 6, height = 4)
F_mean_sq_err_longlong %>%
  # filter(as.numeric(Jitter_radius)<=20) %>%
  mutate(Setup=as.factor(Setup), Jitter_radius=as.factor(as.numeric(Jitter_radius))) %>%
  mutate(Jitter_radius=factor(Jitter_radius,levels=rev(levels(Jitter_radius)))) %>%
  mutate(Setup=recode(Setup, 
                      `1`="Adaptive order", 
                      `3`="Adaptive order", 
                      `2`="Fixed order",
                      `4`="Fixed order",)) %>%
  mutate(Setup=factor(Setup,levels=c('Adaptive order','PPSBM','Fixed order'))) %>%
  ggplot(mapping = aes(x=Jitter_radius, y=F_mse, fill=Setup, color=Setup)) +
  # geom_boxplot() +
  # scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  stat_summary(aes(group=Setup), position = position_dodge(.0),
               geom="point",
               fun = mean) +
  stat_summary(aes(group=Setup),position = position_dodge(.0),
               geom="line",
               fun = "mean") +
  scale_color_manual(values=c('Adaptive order'=scales::hue_pal()(4)[1],'Fixed order'=scales::hue_pal()(4)[2])) +
  xlab("Jitter level") + 
  scale_y_log10()
dev.off()



pdf(file="../Results/Plots/Temp/v_vs_jitter_rad.pdf", width = 6, height = 4)
v_mean_sq_err_longlong %>%
  # filter(as.numeric(Jitter_radius)<=20) %>%
  mutate(Setup=as.factor(Setup), Jitter_radius=as.factor(as.numeric(Jitter_radius))) %>%
  mutate(Jitter_radius=factor(Jitter_radius,levels=rev(levels(Jitter_radius)))) %>%
  mutate(Setup=recode(Setup, 
                      `1`="Adaptive order", 
                      `3`="Adaptive order", 
                      `2`="Fixed order", 
                      `4`="Fixed order")) %>%
  mutate(Setup=factor(Setup,levels=c('Adaptive order','PPSBM','Fixed order'))) %>%
  ggplot(mapping = aes(x=Jitter_radius, y=v_mse, fill=Setup, color=Setup)) +
  # geom_boxplot() +
  # scale_x_discrete(labels = function(x) str_wrap(x, width = 13))
  stat_summary(aes(group=Setup), position = position_dodge(.0),
               geom="point",
               fun = mean) +
  stat_summary(aes(group=Setup),position = position_dodge(.0),
               geom="line",
               fun = "mean") +
  xlab("Jitter level") + 
  scale_y_log10() +
  scale_color_manual(values=c('Adaptive order'=scales::hue_pal()(4)[1],'Fixed order'=scales::hue_pal()(4)[2]))

dev.off()





