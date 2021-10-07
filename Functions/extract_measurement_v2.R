
# Extract certain measurement from a folder_path with same param_setting but various param_value and multiple replicates
# Output: dataframe, cols: param_value | ARI | F_mse | v_mse | 

extract_measurement_v2 = function(folder_path, measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err")){
  
  measurement_df = data.frame()
  
  param_value_vec = list.files(folder_path)
  for (param_value in param_value_vec) {
    file_name_vec = list.files(path = paste0(folder_path,"/",param_value), full.names = T, recursive = TRUE)
    for (file in file_name_vec) {
      load(file)
      meas_value_mat = t(sapply(results, function(one_trial) tryCatch(unlist(one_trial[measurement]), 
                                                                      error=function(x)NA)))
      colnames(meas_value_mat) = measurement
      meas_value_df = as.data.frame(cbind("param_value"=as.numeric(param_value), meas_value_mat))
      measurement_df = dplyr::bind_rows(measurement_df, meas_value_df)
    }
  }
  
  
  return(measurement_df)
}
