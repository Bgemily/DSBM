
# Extract certain measurement from a file_list ---------------------------------


extract_measurement = function(file_list, measurement, obj_name_vec=c("results1")){
  
  measurement_df_list = vector("list",length(obj_name_vec))
  for (f in 1:length(file_list)) {
    file = file_list[[f]]
    load(file)
    
    # if (!is.null(results1$N_node_90[[1]]$network_param$t_vec)) {
    #   t_vec = results1$N_node_90[[1]]$network_param$t_vec
    #   
    # }
    
    for (i in 1:length(obj_name_vec)) {
      obj_name = obj_name_vec[i]
      results1 = get(obj_name)
      tmp_df = c()
      for (j in 1:length(results1)) {
        result_setting_i = results1[[j]]
        # browser()
        tmp = sapply(result_setting_i, function(one_trial) tryCatch(one_trial[[measurement]], 
                                                                    error=function(x)NA))
        tmp_df = cbind(tmp_df, tmp)
      }
      tmp_df = as.data.frame(tmp_df)
      colnames(tmp_df) = names(results1)
      
      #####
      # if (measurement=="F_mean_sq_err") {
      #   tmp_df = tmp_df * (t_vec[2]-t_vec[1])
      # }
      #####
      measurement_df_list[[i]] = rbind(measurement_df_list[[i]], tmp_df)
    }
    
    
    
    
  }
  
  return(measurement_df_list)
}








