
# Extract certain measurement from a folder_path with same param_setting but various param_value and multiple replicates
# Output: dataframe, cols: param_value | ARI | F_mse | v_mse | 

extract_measurement_v2 = function(folder_path, measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err")){
  
  measurement_df = data.frame()
  
  param_value_vec = list.files(folder_path)
  for (param_value in param_value_vec) {
    file_name_vec = list.files(path = paste0(folder_path,"/",param_value), full.names = T, recursive = TRUE)
    for (file in file_name_vec) {
      load(file)
      if (FALSE){
        # Size of meas_value_mat: N_meas*N_trial
        ind = which(measurement=='time_estimation')
        time_est_value_vec = sapply(results, function(one_trial) tryCatch(as.numeric(one_trial[["time_estimation"]],
                                                                                      units = 'secs'), 
                                                                          error=function(x)NA))  
        meas_value_mat = sapply(results[sapply(results,is.list)], function(one_trial) tryCatch(unlist(one_trial[measurement[-ind]]), 
                                                                      error=function(x)NA)) 
        meas_value_mat = rbind(meas_value_mat, time_est_value_vec)
        
      } else{
        # Size of meas_value_mat: N_meas*N_trial
        meas_value_mat = sapply(results[sapply(results,is.list)], function(one_trial) tryCatch(unlist(one_trial[measurement]), 
                                                                      error=function(x)NA))  
      }
      
      if (is.vector(meas_value_mat)) 
        meas_value_mat = matrix(meas_value_mat)
      else
        meas_value_mat = t(meas_value_mat)

      if (ncol(meas_value_mat) == length(measurement)) 
        colnames(meas_value_mat) = measurement
      else if (length(measurement) == 1)
        colnames(meas_value_mat) = rep(measurement,ncol(meas_value_mat))
      
      meas_value_df = as.data.frame(cbind("param_value"=as.numeric(param_value), meas_value_mat))
      measurement_df = dplyr::bind_rows(measurement_df, meas_value_df)
    }
  }
  
  
  return(measurement_df)
}
