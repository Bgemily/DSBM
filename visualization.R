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
library('scales')


# cdf vs pdf vs ppsbm (V!=0) ------------------------------------------------------------
path_vec = rep(0,6)

path_vec[1] = "../Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v19/pr=0.9,n=30,beta=1.3,V=80/"
# path_vec[2] = "../Results/Rdata/SNR_Vnot0_v4/main_v5_pdf_v7_freqtrun7/pr=0.9,n=30,beta=1.3,V=80/"
path_vec[2] = "../Results/Rdata/SNR_Vnot0_v4/main_v5_pdf_v10_freqtrun4/pr=0.9,n=30,beta=1.3,V=80/"
path_vec[3] = "../Results/Rdata/SNR_Vnot0_v4/apply_ppsbm_v4/pr=0.9,n=30,beta=1.3,V=80/"


param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err", 
                                         "N_iteration",
                                         "time_estimation")))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="CDF"), 
                         bind_cols(results_list[[3]],"method"="PPSBM"),
                         # bind_cols(results_list[[4]],"method"="CDF-small-prob-err"),
                         # bind_cols(results_list[[5]],"method"="PDF-small-prob-err"),
                         bind_cols(results_list[[2]],"method"="PDF")
  )
  
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse","V_mse",
                        "N_iteration",
                        "time_per_iter",
                        "time_est")) {
    ### Manipulate column "param_value"
    if (measurement %in% c("1-ARI","f_mse","V_mse")){
      if (param_name == 'n') {
        results_df_tmp = results_df %>%
          filter(param_value <= 90)
      } else{
        results_df_tmp = results_df
      }
      results_df_tmp = results_df_tmp %>%
        mutate(param_value = factor(param_value, levels = sort(unique(param_value),
                                                               decreasing = param_name=="V")) )
    } else{
      results_df_tmp = results_df 
    }
    results_df_tmp = results_df_tmp %>% 
      mutate(time_per_iter=time_estimation/N_iteration)
    pdf(file=paste0("../Results/Plots/Temp/Cdf_vs_pdf_vs_ppsbm_Vnot0/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'vs', '_', 
                    switch(param_name,"beta"='beta','n'='p','V'='V'),
                    ".pdf"), 
        width = 4, height = 2)
    g = results_df_tmp %>% 
      mutate(time_per_iter=time_estimation/N_iteration,
             method=as_factor(method)) %>%
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "time_est"= time_estimation,
                          "N_iteration" = N_iteration,
                          "time_per_iter" = time_per_iter,
                          "V_mse" = v_mean_sq_err), 
                 # linetype=method,
                 color=method)) +
      stat_summary(aes(group=method),position = position_dodge(.0),
                   geom="line",
                   size=0.7,
                   alpha=0.8,
                   fun = "mean") +
      # stat_summary(data = results_df_tmp %>%
      #                filter(param_value==ifelse(is.factor(results_df_tmp$param_value),
      #                                           levels(results_df_tmp$param_value)[1],
      #                                           sort(results_df_tmp$param_value, 
      #                                                decreasing = param_name=="V")[1])),
      #              aes(group=method), position = position_dodge(.0),
      #              geom="point",
      #              fun = "mean",
      #              alpha=0.7) +
      theme_bw()+
      theme(legend.position = "right", 
            legend.box.background = element_rect(color="gray")) +
      guides(color=guide_legend(title="Method")) +
      scale_color_manual(breaks=c("CDF","PDF",'PPSBM'),
                         values=hue_pal()(3), 
                         labels=c("SibSBM-C","SibSBM-P",'PPSBM'),
                         drop=FALSE) +
      # guides(color=guide_legend(nrow=2,byrow=TRUE),
      #        linetype=guide_legend(nrow=2,byrow=TRUE)) +
      coord_cartesian(ylim = switch(measurement,
                                    "N_iteration" = c(0,11),
                                    "time_est" = c(0,400),
                                    "time_per_iter" = c(0,20),
                                    "V_mse" = c(0,150),
                                    "f_mse" = c(0,0.04),
                                    "1-ARI" = c(0,1) )) +
      ylab(switch(measurement, 'time_est'='Total cmpt time (s)',
                  "time_per_iter"='Cmpt time per iter (s)',
                  "N_iteration"='Number of iterations',
                  measurement) ) +
      xlab(switch(param_name, 'n'="p",
                  'beta'="beta",
                  'V'="V"))
    
    print(g)
    dev.off()
  }
  
  
}




# cdf vs pdf vs ppsbm (V==0) ------------------------------------------------------------
path_vec = rep(0,6)

path_vec[1] = "../Results/Rdata/SNR_Vis0_v4/main_v5_cdf_v3/pr=0.9,n=30,beta=1.3,V=0/"
path_vec[2] = "../Results/Rdata/SNR_Vis0_v4/main_v5_pdf_v6_freqtrun7/pr=0.9,n=30,beta=1.3,V=0/"
path_vec[3] = "../Results/Rdata/SNR_Vis0_v4/apply_ppsbm_v4/pr=0.9,n=30,beta=1.3,V=0/"


param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err", 
                                         "N_iteration",
                                         "time_estimation")))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="CDF"), 
                         bind_cols(results_list[[3]],"method"="PPSBM"),
                         bind_cols(results_list[[2]],"method"="PDF")
  )
  
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse",
                        "N_iteration",
                        "time_per_iter",
                        "time_est"
                        )) {
    ### Manipulate column "param_value"
    if (measurement %in% c("1-ARI","f_mse","V_mse")){
      if (param_name == 'n') {
        results_df_tmp = results_df %>%
          filter(param_value <= 90)
      } else{
        results_df_tmp = results_df
      }
      results_df_tmp = results_df_tmp %>%
        mutate(param_value = factor(param_value, levels = sort(unique(param_value),
                                                               decreasing = param_name=="V")) )
    } else{
      results_df_tmp = results_df
    }
    
    pdf(file=paste0("../Results/Plots/Temp/Cdf_vs_pdf_vs_ppsbm_Vis0/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'vs', '_', 
                    switch(param_name,"beta"='beta','n'='p','V'='V'),
                    ".pdf"), 
        width = 4, height = 2)
    g = results_df_tmp %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "time_est"= time_estimation,
                          "N_iteration" = N_iteration,
                          "time_per_iter" = time_estimation/N_iteration,
                          "V_mse" = v_mean_sq_err), 
                 # linetype=method,
                 color=method)) +
      stat_summary(aes(group=method),position = position_dodge(.0),
                   geom="line",
                   alpha=0.8,
                   fun = "mean") +
      stat_summary(data = results_df_tmp %>%
                     filter(param_value==ifelse(is.factor(results_df_tmp$param_value),
                                                levels(results_df_tmp$param_value)[1],
                                                sort(results_df_tmp$param_value, 
                                                     decreasing = param_name=="V")[1])),
                   aes(group=method), position = position_dodge(.0),
                   geom="point",
                   fun = "mean",
                   alpha=0.7) +
      theme_bw()+
      theme(legend.position = "right", 
            legend.box.background = element_rect(color="gray")) +
      guides(color=guide_legend(title="Method"))+
      # guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      coord_cartesian(ylim = switch(measurement,
                                    # "time_est" = c(0,350),
                                    "V_mse" = c(0,150),
                                    "f_mse" = c(0,0.016),
                                    "1-ARI" = c(0,1)
                                    )) +
      ylab(switch(measurement, 'time_est'='Total cmpt time (s)',
                  "time_per_iter"='Cmpt time per iter (s)',
                  "N_iteration"='Number of iterations',
                  measurement) ) +
      xlab(switch(param_name, 'n'='N_node',param_name))
    if (measurement=='time_est') {
      g = g +
        scale_x_continuous(limits = c(30,150), breaks=c(30,60,90,120,150))+
        scale_y_continuous(limits = c(0.2,250), breaks = c(0.2,100,200,250))
        # scale_y_continuous(trans = 'log10')
    }
    print(g)
    dev.off()
  }
  
}



# ICL vs N_clus (separate PDF and CDF in different plots) -----------------------------

path_vec = rep(0,2)

for (method in c("CDF","PDF")) {
  if(method=='CDF'){
    freq_trun_vec = "NA"
    path_vec[1] = "../Results/Rdata/ICL_v1/main_v5_cdf_v2/pr=0.9,n=30,beta=1.9,V=80/"
  } else{
    freq_trun_vec = 1:6
    path_vec[1] = "../Results/Rdata/ICL_v1/main_v5_pdf_v4/pr=0.9,n=30,beta=1.9,V=80//"
  }
  
  param_name_vec = list.files(path_vec[1])
  
  ### For each parameter (n/beta/V), extract results and visualize results
  for (param_name in param_name_vec) {
    
    ### Extract results for n/beta/V. Output: freq_trun | ICL | log-lik | penalty
    func_tmp = function(folder_path, param_name) {
      extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name),
                             measurement = c("ICL_mat","N_clus_est"))
    }
    results_list = lapply(path_vec, func_tmp, param_name=param_name)
    results_df = bind_rows(bind_cols(results_list[[1]][,-ncol(results_list[[1]])])
    )
    results_df = colMeans(results_df[,-1])
    N_clus = 5
    results_df = as.tibble(matrix(results_df, ncol=N_clus))
    colnames(results_df) = paste0("N_clus",1:5)
    results_df = bind_cols(results_df, freq_trun=freq_trun_vec) %>% mutate(freq_trun = as_factor(freq_trun))
    results_df = results_df %>% 
      # filter(param_value==1.9) %>%
      # mutate(param_value = as_factor(param_value), freq_trun = as_factor(freq_trun)) %>%
      pivot_longer(cols = starts_with("N_clus"), names_to = "N_clus_ICL", values_to = "ICL")
    
    
    ### Plot ICL vs N_clus vs freq_trun
    for (measurement in c("ICL")) {
      pdf(file=paste0("../Results/Plots/Temp/ICL_vs_Nclus/",
                      if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                      '_', 'N_clus', '_', method, ".pdf"),
          width = 4, height = 2)
      g = results_df %>%
        ggplot(aes(x=switch(measurement,
                            "ICL" = N_clus_ICL,
                            "compl_log_lik" = N_clus_compl_log_lik,
                            "log_lik" = N_clus_log_lik,
                            "penalty_2" = N_clus_penalty_2,
                            "penalty" = N_clus_penalty),
                   y=switch(measurement,
                            "ICL" = ICL,
                            "compl_log_lik" = compl_log_lik,
                            "log_lik" = log_lik,
                            "penalty_2" = penalty_2,
                            "penalty" = penalty),
                   size=freq_trun,
                   alpha=freq_trun,
                   color=freq_trun)) +
        stat_summary(aes(group=freq_trun),position = position_dodge(.0),
                     geom="line",
                     # alpha=0.7,
                     fun = "mean") +
        stat_summary(aes(group=freq_trun),position = position_dodge(.0),
                     geom="point",
                     # alpha=0.7,
                     size=1,
                     fun = "mean") +
        scale_color_manual(values = c(hue_pal()(3)[1],
                                      brewer_pal(palette = "Dark2")(8)[2:3],
                                      hue_pal()(3)[2],
                                      brewer_pal(palette = "Dark2")(8)[c(4,6)]
                                      ),
                           breaks = c('NA',2:6)) +
        scale_size_manual(values = c(0.7, 
                                      rep(0.5,2),
                                      0.7,
                                      rep(0.5,2)),
                            breaks = c('NA',2:6)) +
        scale_alpha_manual(values = c(0.9, 
                                     rep(0.4,2),
                                     0.9,
                                     rep(0.4,2)),
                          breaks = c('NA',2:6)) +
        theme(legend.position = "right") +
        guides(color=guide_legend(title=element_blank()),
               alpha=guide_legend(title=element_blank()),
               size=guide_legend(title=element_blank())) +
        coord_cartesian(ylim = switch(method,"PDF"=c(-2000,-1700))) +
        ylab(measurement) +
        scale_x_discrete(labels=1:5) +
        xlab("Number of clusters")
      
      print(g)
      dev.off()
    }
  }
  
}


# ICL vs N_clus (put PDF and CDF on same plot) ----------------------------

path_vec = rep(0,2)
path_vec[1] = "../Results/Rdata/ICL_v1/main_v5_cdf_v2/pr=0.9,n=30,beta=1.6,V=80/"
path_vec[2] = "../Results/Rdata/ICL_v1/main_v5_pdf_v4/pr=0.9,n=30,beta=1.7,V=80//"
param_name_vec = list.files(path_vec[1])

freq_trun_vec_CDF = c("NA")
freq_trun_vec_PDF = 1:6

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: freq_trun | ICL | log-lik | penalty
  func_tmp = function(folder_path, param_name) {
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name),
                           measurement = c("ICL_mat","N_clus_est"))
  }
  results_list = lapply(path_vec, func_tmp, param_name=param_name)
  
  results_df_CDF = bind_rows(bind_cols(results_list[[1]][,-ncol(results_list[[1]])])
  )
  results_df_CDF = colMeans(results_df_CDF[,-1])
  N_clus = 5
  results_df_CDF = as.tibble(matrix(results_df_CDF, ncol=N_clus))
  colnames(results_df_CDF) = paste0("N_clus",1:5)
  results_df_CDF = bind_cols(results_df_CDF, freq_trun=freq_trun_vec_CDF) %>% mutate(freq_trun = as_factor(freq_trun))
  results_df_CDF = results_df_CDF %>% 
    pivot_longer(cols = starts_with("N_clus"), names_to = "N_clus_ICL", values_to = "ICL")
  
  results_df_PDF = bind_rows(bind_cols(results_list[[2]][,-ncol(results_list[[2]])])
  )
  results_df_PDF = colMeans(results_df_PDF[,-1])
  N_clus = 5
  results_df_PDF = as.tibble(matrix(results_df_PDF, ncol=N_clus))
  colnames(results_df_PDF) = paste0("N_clus",1:5)
  results_df_PDF = bind_cols(results_df_PDF, freq_trun=freq_trun_vec_PDF) %>% mutate(freq_trun = as_factor(freq_trun))
  results_df_PDF = results_df_PDF %>% 
    pivot_longer(cols = starts_with("N_clus"), names_to = "N_clus_ICL", values_to = "ICL")
  
  results_df = bind_rows(results_df_CDF, results_df_PDF)
  
  ### Plot ICL vs N_clus vs freq_trun
  for (measurement in c("ICL")) {
    pdf(file=paste0("../Results/Plots/Temp/ICL_vs_Nclus/",
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'N_clus', ".pdf"),
        width = 3.8, height = 3.5)
    g = results_df %>%
      filter(freq_trun %in% c("NA",2,3,4)) %>%
      ggplot(aes(x=switch(measurement,
                          "ICL" = N_clus_ICL,
                          "compl_log_lik" = N_clus_compl_log_lik,
                          "log_lik" = N_clus_log_lik,
                          "penalty_2" = N_clus_penalty_2,
                          "penalty" = N_clus_penalty),
                 y=switch(measurement,
                          "ICL" = ICL,
                          "compl_log_lik" = compl_log_lik,
                          "log_lik" = log_lik,
                          "penalty_2" = penalty_2,
                          "penalty" = penalty),
                 size=freq_trun,
                 alpha=freq_trun,
                 color=freq_trun)) +
      stat_summary(aes(group=freq_trun),position = position_dodge(.0),
                   geom="line", size=0.7,
                   # alpha=0.7,
                   fun = "mean") +
      # stat_summary(aes(group=freq_trun),position = position_dodge(.0),
      #              geom="point",
      #              # alpha=0.7,
      #              size=1,
      #              fun = "mean") +
      scale_color_manual(values = c(hue_pal()(3)[1],
                                    brewer_pal(palette = "Dark2")(8)[c(3)],
                                    hue_pal()(3)[2],
                                    brewer_pal(palette = "Dark2")(8)[c(4)]),
                         labels = c("SidSBM-C", paste0("SidSBM-P","(    = ",2:4,")")),
                         breaks = c('NA',2:4)) +
      scale_size_manual(values = c(0.7, 
                                   0.5,
                                   0.7,
                                   0.5),
                        labels = c("SidSBM-C", paste0("SidSBM-P","(    = ",2:4,")")),
                        breaks = c('NA',2:4)) +
      scale_alpha_manual(values = c(0.9, 
                                    0.4,
                                    0.9,
                                    0.4),
                         labels = c("SidSBM-C", paste0("SidSBM-P","(    = ",2:4,")")),
                         breaks = c('NA',2:4)) +
      theme_bw()+
      theme(legend.position = c(0.99,0.99),
            legend.justification = c("right","top"),
            legend.background = element_rect(fill='transparent'),
            legend.box.background = element_rect(color="gray"),
            legend.key = element_rect(size=1, fill = NA, color = NA),
            legend.key.height = unit(.4, 'cm'),
            # legend.key.width = unit(.5, 'cm'),
            legend.title = element_text(size=9.5, 
                                        margin=margin(b=-3,unit='pt')),
            legend.margin = margin(3,3,3,3),
            legend.text = element_text(size = 8.5 )) +
      guides(color=guide_legend(title="Method", nrow=2, byrow=F),
             alpha=guide_legend(title="Method", nrow=2, byrow=F),
             size=guide_legend(title="Method", nrow=2, byrow=F)) +
      coord_cartesian(ylim = c(-2100,-1250)) +
      ylab(measurement) + 
      scale_x_discrete(labels=1:5) +
      xlab("Number of clusters")
    
    print(g)
    dev.off()
  }
}



# ARI/F_mse/v_mse vs iteration number -------------------------------------
path_vec = rep(0,6)

path_vec[1] = "../Results/Rdata/Initialization_v1/main_v5_cdf_v2/Our_init/pr=0.9,n=30,beta=1.9,V=80/"
path_vec[2] = "../Results/Rdata/Initialization_v1/main_v5_cdf_v2/N_restart_v2/1/pr=0.9,n=30,beta=1.9,V=80/"
path_vec[3] = "../Results/Rdata/Initialization_v1/main_v5_cdf_v2/N_restart_v2/10/pr=0.9,n=30,beta=1.9,V=80/"
path_vec[4] = "../Results/Rdata/Initialization_v1/main_v5_cdf_v2/N_restart_v2/2/pr=0.9,n=30,beta=1.9,V=80/"
path_vec[5] = "../Results/Rdata/Initialization_v1/main_v5_cdf_v2/N_restart_v2/3/pr=0.9,n=30,beta=1.9,V=80/"

param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_history", "ARI_mean",
                                         "F_mean_sq_err_history", 
                                         'time_init',
                                         "v_mean_sq_err_history")))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="Proposal"),
                         # bind_cols(results_list[[3]],"method"="10 random"),
                         bind_cols(results_list[[4]],"method"="2 random"),
                         bind_cols(results_list[[5]],"method"="3 random"),
                         bind_cols(results_list[[2]],"method"="1 random"))
  
  results_df = results_df %>% 
    mutate(param_value = as_factor(param_value)) %>%
    pivot_longer(cols = starts_with("ARI_history"), names_to = "Iter_ARI", values_to = "ARI") %>% 
    pivot_longer(cols = starts_with("F_mean_sq_err"), names_to = "Iter_F_mean_sq_err", values_to = "F_mean_sq_err") %>% 
    pivot_longer(cols = starts_with("v_mean_sq_err"), names_to = "Iter_v_mean_sq_err", values_to = "v_mean_sq_err")
  
  results_df = results_df %>% 
    mutate(method=fct_relevel(method,"10 random", "20 random", after=Inf),
           method=fct_relevel(method,"Proposal", after=0),
           Iter_ARI=fct_relevel(Iter_ARI,"ARI_history10","ARI_history11", after = Inf),
           Iter_F_mean_sq_err=fct_relevel(Iter_F_mean_sq_err,"F_mean_sq_err_history10","F_mean_sq_err_history11", after = Inf),
           Iter_v_mean_sq_err=fct_relevel(Iter_v_mean_sq_err,"v_mean_sq_err_history10","v_mean_sq_err_history11", after = Inf))
  
  ### Plot ICL vs N_clus
  for (measurement in c("F_mean_sq_err",'v_mean_sq_err',"1-ARI")) {
    pdf(file=paste0("../Results/Plots/Temp/Initialization/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement), 
                    '_', 'Iteration', ".pdf"), 
        width = 3.5, height = 2)
    g = results_df %>% 
      filter(Iter_ARI!="ARI_history11",
             Iter_F_mean_sq_err!="F_mean_sq_err_history11",
             Iter_v_mean_sq_err!="v_mean_sq_err11") %>%
      ggplot(aes(x=switch(measurement,
                          "1-ARI" = Iter_ARI,
                          "F_mean_sq_err" = Iter_F_mean_sq_err, 
                          "v_mean_sq_err" = Iter_v_mean_sq_err),
                 y=switch(measurement,
                          "1-ARI" = 1-ARI,
                          "F_mean_sq_err" = F_mean_sq_err,
                          "v_mean_sq_err" = v_mean_sq_err),
                 linetype=method,
                 color=method)) +
      stat_summary(aes(group=method),position = position_dodge(.0),
                   geom="line",alpha=0.9, size=0.7,
                   fun = switch(measurement,
                                "F_mean_sq_err"="median",
                                "v_mean_sq_err"="median",
                                "mean")) +
      theme_bw()+
      theme(legend.position = switch(measurement, "1-ARI"=c(0.995,0.99), "none"),
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill='transparent'),
            legend.box.background = element_rect(color="gray"),
            legend.key = element_rect(size=1, fill = NA, color = NA),
            legend.key.height = unit(.4, 'cm'),
            legend.key.width = unit(.8, 'cm'),
            legend.title = element_text(size=9.5, 
                                        margin=margin(b=-3,unit='pt')),
            legend.margin = margin(3,3,3,3),
            legend.text = element_text(size = 9.5 )) +
      scale_color_discrete(breaks = c("Proposal",
                                      "1 random",
                                      "2 random",
                                      "3 random",
                                      "5 random",
                                      "10 random",
                                      "20 random"),
                         labels = c("Proposal",
                                    "Random (1 restart)",
                                    "Random (2 restarts)",
                                    "Random (3 restarts)",
                                    "Random (5 restarts)",
                                    "Random (10 restarts)",
                                    "Random (20 restarts)")) +
      scale_linetype_manual(breaks = c("Proposal",
                                       "1 random",
                                       "2 random",
                                       "3 random",
                                       "5 random",
                                       "10 random",
                                       "20 random"),
                            labels = c("Proposal",
                                       "Random (1 restart)",
                                       "Random (2 restarts)",
                                       "Random (3 restarts)",
                                       "Random (5 restarts)",
                                       "Random (10 restarts)",
                                       "Random (20 restarts)"),
                            values = c("solid",rep('dashed',6))) +
      guides(color=guide_legend(title="Initialization Scheme"), 
             linetype=guide_legend(title="Initialization Scheme")) +
      coord_cartesian(ylim=switch(measurement,
                                  "1-ARI"=c(0,1),
                                  "F_mean_sq_err" = c(0,8),
                                  "v_mean_sq_err" = c(0,100)
                                  )) +
      ylab(switch(measurement,
                  "F_mean_sq_err"="F_mse",
                  measurement)) +
      scale_x_discrete(labels=c(0:10)) +
      xlab("Iteration")
    
    print(g)
    dev.off()
  }
  
}



# ARI vs Iteration for n=30,66,150 ----------------------------------------

load("/Users/ztzhang/Documents/Academic/SC/graphon/Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v19/pr=0.9,n=30,beta=1.3,V=80/n/30/N_trial10_20220113_094142.Rdata")
lapply(results,'[[','ARI_history')->ARI_history_list
lapply(results,'[[','F_mean_sq_err_history')->F_mean_sq_err_history_list

load("/Users/ztzhang/Documents/Academic/SC/graphon/Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v19/pr=0.9,n=30,beta=1.3,V=80/n/66/N_trial10_20220113_094434.Rdata")
lapply(results,'[[','F_mean_sq_err_history')->F_mean_sq_err_history_list_2

load("/Users/ztzhang/Documents/Academic/SC/graphon/Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v19/pr=0.9,n=30,beta=1.3,V=80/n/150/N_trial10_20220113_095150.Rdata")
lapply(results,'[[','F_mean_sq_err_history')->F_mean_sq_err_history_list_3

plot(0:10,rep(0,11),type='l',ylim=c(0,1),col=0,xlab = 'Iteration',ylab="1-ARI",
     main = 'Black/Red/Green: 30/66/150 nodes')
lines(1-ARI_history_list[[3]])
lines(1-ARI_history_list_2[[4]],col=2)
lines(1-ARI_history_list_3[[2]],col=3)

# Archive -----------------------------------------------------------------


# cdf vs pdf  ------------------------------------------------------------
path_vec = rep(0,6)

path_vec[1] = "../Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v11_freqtrun7/pr=0.9,n=30,beta=1.3,V=80/"
path_vec[2] = "../Results/Rdata/SNR_Vnot0_v4/main_v5_pdf_v3_freqtrun7/pr=0.9,n=30,beta=1.3,V=80/"


param_name_vec = list.files(path_vec[1])

### For each parameter (n/beta/V), extract results and visualize results
for (param_name in param_name_vec) {
  
  ### Extract results for n/beta/V. Output: param_value (n/beta/V's value) | ARI | F_mse | V_mse | method
  results_list = lapply(path_vec, function(folder_path)
    extract_measurement_v2(folder_path = paste0(folder_path,"/",param_name), 
                           measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err", 
                                         "N_iteration",
                                         "time_estimation")))
  results_df = bind_rows(bind_cols(results_list[[1]],"method"="CDF"), 
                         # bind_cols(results_list[[3]],"method"="PPSBM"),
                         # bind_cols(results_list[[4]],"method"="PDF+N_basis_11"),
                         # bind_cols(results_list[[5]],"method"="PDF+N_basis_19"),
                         bind_cols(results_list[[2]],"method"="PDF")
                         )
  
  ### Manipulate column "param_value"
  results_df = results_df %>%
    mutate(param_value = factor(param_value, levels = sort(unique(param_value),
                                                           decreasing = param_name=="V")) )
  
  ### Plot ARI/F_mse vs n/beta/V
  for (measurement in c("1-ARI","f_mse","V_mse",
                        "N_iteration",
                        "time_per_iter",
                        "time_est")) {
    pdf(file=paste0("../Results/Plots/Temp/", 
                    if_else(measurement=="1-ARI", true = "ARI", false = measurement),
                    '_', 'vs', '_', 
                    switch(param_name,"beta"='beta','n'='p','V'='V'),
                    ".pdf"), 
        width = 5, height = 5)
    g = results_df %>% 
      ggplot(aes(x=param_value, 
                 y=switch(measurement,
                          "1-ARI" = 1-ARI_mean,
                          "f_mse" = F_mean_sq_err,
                          "time_est"= time_estimation,
                          "N_iteration" = N_iteration,
                          "time_per_iter" = time_estimation/N_iteration,
                          "V_mse" = v_mean_sq_err), 
                 linetype=method,
                 color=method)) +
      stat_summary(aes(group=method), position = position_dodge(.0),
                   # geom="pointrange", 
                   alpha=0.7,
                   # fun.min = function(x)quantile(x,0.25),
                   # fun.max = function(x)quantile(x,0.75),
                   # fun.min = function(x)mean(x)-sd(x),
                   # fun.max = function(x)mean(x)+sd(x),
                   geom="point",
                   fun = "mean") +
      stat_summary(aes(group=method),position = position_dodge(.0),
                   geom="line",
                   alpha=0.8,
                   fun = "mean") +
      theme(legend.position = "bottom") +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
      scale_y_continuous(limits = switch(measurement,
                                         # "time_est" = c(0,350),
                                         "V_mse" = c(0,150),
                                         "f_mse" = c(0,0.025),
                                         "1-ARI" = c(0,1) )) +
      ylab(measurement) +
      xlab(switch(param_name, 'n'='N_node',param_name))
    
    print(g)
    dev.off()
  }
  
}



# pdf performance vs N_node (V!=0) ------------------------------------------------------------
path_vec = rep(0,5)

path_vec[1] = "../Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v2_multi_Nclus/pr=0.9,n=60,beta=1.9,V=80,N_grid=20/"

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





