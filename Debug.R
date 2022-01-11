
# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# -------------------------------------------------------------------------
load("../Results/Rdata/SNR_Vnot0_v4/main_v5_cdf_v11_freqtrun7/pr=0.9,n=30,beta=1.3,V=80/beta/1.3/N_trial10_20211206_231811.Rdata")

network_params = do.call(generate_network2_v3, args=results[[1]]$network_param)
network_params$time_shift_list

plot(network_params$time_shift_list[[1]], y=results[[1]]$v_vec_list_est[[1]])
abline(a=0,b=1,col=2)
