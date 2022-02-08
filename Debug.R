
# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# -------------------------------------------------------------------------

res = 3
for(res in 1:10){
  v_est = results[[res]]$v_vec_list_est[[1]]
  generated_network = do.call(generate_network2_v3, results[[res]]$network_param)
  v_true = generated_network$time_shift_list[[1]]
  
  plot(v_true,v_est,main=round(results[[res]]$v_mean_sq_err,2)); abline(a=0,b=1,col=2)
  results[[res]]$v_mean_sq_err
  Sys.sleep(1)
}



