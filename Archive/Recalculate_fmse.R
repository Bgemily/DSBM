recalculate_fmse= function(file_list){
  
  F_mean_sq_err_longlong = data.frame()
  
  for (file in file_list) {
    load(file)
    results_vec = c("results1", "results2")
    if (!is.null(results1$N_node_90[[1]]$t_vec)) {
      t_vec = results1$N_node_90[[1]]$t_vec
      
      F_mean_sq_err_list = extract_measurement(file_list = list(file), 
                                               measurement = "F_mean_sq_err", 
                                               obj_name = results_vec)
      
      for (i in 1:length(results_vec)) {
        obj_name = results_vec[i]
        F_mean_sq_err_longlong = data.frame()
        for (j in 1:length(path_vec)) {
          file_list = file_list_mat[[j]]
          
          F_mean_sq_err = F_mean_sq_err_list[[j]][[i]]
          F_mean_sq_err = rename_with(F_mean_sq_err, function(s)as.numeric(gsub("[^0-9.]", "", s)))
          F_mean_sq_err_long = F_mean_sq_err %>% 
            pivot_longer(cols = everything(), 
                         names_to = "Jitter_radius", 
                         values_to = "F_mse") %>%
            mutate(Setup = j)
          F_mean_sq_err_longlong = rbind(F_mean_sq_err_longlong, F_mean_sq_err_long)
          
        }
      }
        
    }
  }
  
  
  
  
  return(F_mean_sq_err_longlong)
}