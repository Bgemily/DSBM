
# Implement the algorithm adapted from Bigot and Gendre (2013)

##### Algorithm #####
# Input: Edge time matrix list, clusters list, step size, t_vec
# Output: n0_vec


est_n0_vec_v3 = function(edge_time_mat_list, clusters_list, t_vec=seq(0,50,0.05), 
                        step_size=0.02, max_iter=100, epsilon=0.001, 
                        n0_vec_list_true){
  
  t_unit = t_vec[2] - t_vec[1]
  N_subj = length(edge_time_mat_list)
  N_clus = length(clusters_list[[1]])
  N_node_list = lapply(edge_time_mat_list, nrow)
  n0_vec_list = lapply(N_node_list, numeric)
  membership_list = lapply(clusters_list, clus2mem)
  
  for (m in 1:N_subj) {
    diag(edge_time_mat_list[[m]]) = NA
  }
  
  ### Initialize n0_vec by earliest edge time
  earliest_edge_time_list = list()
  for (edge_time_mat in edge_time_mat_list) {
    earliest_edge_time_vec = apply(edge_time_mat, 1, function(row)min(row[which(row>1)]))
    earliest_edge_time_list = c(earliest_edge_time_list, list(earliest_edge_time_vec))
  }
  n0_vec_list = lapply(earliest_edge_time_list, function(v_vec)(v_vec-min(v_vec))/t_unit)
  n0_vec_list = lapply(n0_vec_list, round)
  
  ### Initialize n0_vec by zero
  n0_vec_list = lapply(N_node_list, numeric)
  
  ### Obtain the order of time shifts
  n0_order_list = lapply(n0_vec_list, order)
  
  
  ### Gradient descent
  n_iter = 0
  converge = FALSE
  n0_vec_list_update = n0_vec_list_current = n0_vec_list
  loss_history = c()
  accuracy_history = c()
  while (!converge && n_iter<= max_iter) {
    ### Compute hat_F_array, i.e. hat_F_qk's
    shifted_edge_time_mat_list = list()
    for (m in 1:N_subj) {
      shifted_edge_time_mat = edge_time_mat_list[[m]] - n0_vec2mat(n0_vec_list_current[[m]])*t_unit
      diag(shifted_edge_time_mat) = NA
      shifted_edge_time_mat_list[[m]] = shifted_edge_time_mat
    }
    hat_F_array = array(dim=c(N_clus, N_clus, length(t_vec)))
    for (q in 1:N_clus) {
      for (k in 1:N_clus) {
        tmp = lapply(1:N_subj, function(m) shifted_edge_time_mat_list[[m]][clusters_list[[m]][[q]], clusters_list[[m]][[k]]])
        tmp = unlist(tmp)
        hat_F_qk_vec = ecdf(tmp)(t_vec)
        hat_F_array[q,k,] = hat_F_qk_vec
      }
    }
    
    ### Update n0_vec_list_update
    gradient_vec_list = lapply(N_node_list, numeric)
    for (m in 1:N_subj) {
      for (order_ind in N_node_list[[m]]:1) {
        ### V1
        ### V2: coordinate descent
        # n0_vec_list_current = n0_vec_list_update
        ### ### ###
        
        
        i = n0_order_list[[m]][order_ind]

        ### Get min_n0 and max_n0
        ### V1: force the order to be unchanged
        # if(order_ind>=2 & order_ind<N_node_list[[m]]){
        #   pre_ind = n0_order_list[[m]][order_ind-1]
        #   min_n0 = n0_vec_list_current[[m]][pre_ind]
        #   post_ind = n0_order_list[[m]][order_ind+1]
        #   max_n0 = n0_vec_list_current[[m]][post_ind]
        # }
        # else if(order_ind==1){
        #   min_n0 = -max(t_vec)/t_unit
        #   post_ind = n0_order_list[[m]][order_ind+1]
        #   max_n0 = n0_vec_list_current[[m]][post_ind]
        # }
        # else{
        #   pre_ind = n0_order_list[[m]][order_ind-1]
        #   min_n0 = n0_vec_list_current[[m]][pre_ind]
        #   max_n0 = max(t_vec)/t_unit
        # }
        ### V2: order can change
        max_n0 = (1/3)*max(t_vec)/t_unit
        min_n0 = -(1/3)*max(t_vec)/t_unit
        ### ### ###
        
        
        ### Compute gradient
        q = membership_list[[m]][i] #z_mi
        A = matrix(nrow=N_clus, ncol=length(t_vec)-1)
        B = matrix(nrow=N_clus, ncol=length(t_vec)-1)
        C = vector(length=N_clus)
        for (k in 1:N_clus) {
          ### Define frequency vector ###
          l_vec = 1:(length(t_vec)-1)
          l_vec = c(tail(l_vec, length(l_vec)/2)-length(t_vec),
                    head(l_vec, length(l_vec)/2))
          ### ### ###
          
          ### Compute A[k,] ###
          J_mki = which(membership_list[[m]]==k & n0_vec_list_current[[m]]<=n0_vec_list_current[[m]][i])
          J_mki = setdiff(J_mki, i)
          if (length(J_mki) > 0){
            DdotQ_mki_vec = -length(J_mki) * ecdf(edge_time_mat_list[[m]][i,J_mki])(t_vec)
          }
          else{
            DdotQ_mki_vec = rep(0, length(t_vec))
          }
          
          fft_DdotQ_mki_vec = fft(DdotQ_mki_vec)[-1]
          fft_DdotQ_mki_vec = c(tail(fft_DdotQ_mki_vec, length(fft_DdotQ_mki_vec)/2), 
                                head(fft_DdotQ_mki_vec, length(fft_DdotQ_mki_vec)/2))
          A[k,] = 1i*2*pi*(l_vec/length(t_vec)) * 
                  exp(1i*2*pi*n0_vec_list_current[[m]][i]*(l_vec/length(t_vec))) * 
                  ( -fft_DdotQ_mki_vec - tail(DdotQ_mki_vec,1)/(1-exp(-1i*2*pi*(l_vec/length(t_vec)))) )
          ### ###
          
          ### Compute B[k,] ###
          ### V1
          ### BUG: hat_F_qk_vec should NOT depend on m!
          # shifted_edge_time_mat = edge_time_mat_list[[m]] - n0_vec2mat(n0_vec_list_current[[m]])*t_unit
          # diag(shifted_edge_time_mat) = NA
          # hat_F_qk_vec = ecdf(shifted_edge_time_mat[clusters_list[[m]][[q]], clusters_list[[m]][[k]]])(t_vec)
          ### V2
          hat_F_qk_vec = hat_F_array[q,k,]
          ### ### ###
          
          ### Take fft, and remove the Fourier coef corresponding to l=0 
          fft_hat_F_qk_vec = fft(hat_F_qk_vec)[-1] 
          ### Reorder the Fourier coef s.t. l=-(N-1)/2, \cdots,(N-1)/2, l\neq 0.
          fft_hat_F_qk_vec = c(tail(fft_hat_F_qk_vec, length(fft_hat_F_qk_vec)/2), 
                               head(fft_hat_F_qk_vec, length(fft_hat_F_qk_vec)/2))
          B[k,] = fft_hat_F_qk_vec + 1/(1-exp(-1i*2*pi*(l_vec/length(t_vec))))
          ### ###
          
          ### Compute C[k] ###
          C[k] = fft(hat_F_qk_vec)[1]*tail(DdotQ_mki_vec,1) -
            fft(DdotQ_mki_vec)[1] -
            n0_vec_list_current[[m]][i]*tail(DdotQ_mki_vec,1)
          ### ### ###
          
        } 
        
        gradient = 4 * ( -Re(sum(A*Conj(B))) + sum(C) )
        gradient = Re(gradient) # Get rid of "+0i"
        gradient_vec_list[[m]][i] = gradient
        
        ### Update n0_vec_list
        n0_vec_list_update[[m]][i] = n0_vec_list_current[[m]][i] - step_size*gradient  
        n0_vec_list_update[[m]][i] = round(n0_vec_list_update[[m]][i])
        if(n0_vec_list_update[[m]][i] < min_n0){
          n0_vec_list_update[[m]][i] = min_n0
        }
        else if (n0_vec_list_update[[m]][i] > max_n0){
          n0_vec_list_update[[m]][i] = max_n0
        }
        
        
      }
    }
    
    n0_vec_list_update = lapply(n0_vec_list_update, function(n0_vec) n0_vec-min(n0_vec))
    
    
    ### Evaluate loss function (up to a constant) (this step is time-consuming)
    loss = eval_loss(edge_time_mat_list = edge_time_mat_list, n0_vec_list = n0_vec_list_update,
                     clusters_list = clusters_list, t_vec = t_vec)$loss
    loss_history = c(loss_history, loss)
    
    
    ### Evaluate accuracy of n0_vec_list
    ### V1: Accuracy is correlation
    # accuracy = cor(unlist(n0_vec_list_update), unlist(n0_vec_list_true))
    ### V2: Accuracy is mean squared error
    accuracy = mean( (unlist(n0_vec_list_update)*t_unit - unlist(n0_vec_list_true)*t_unit)^2 )
    ### ### ###
    accuracy_history = c(accuracy_history, accuracy)
    
    ### Evaluate convergence criterion
    converge = sqrt(sum((unlist(n0_vec_list_update)-unlist(n0_vec_list_current))^2))/
      sqrt(sum((unlist(n0_vec_list_current)+.Machine$double.eps)^2)) < epsilon
    
    ### Debug
    # converge = FALSE
    
    n_iter = n_iter + 1
    n0_vec_list_current = n0_vec_list_update
  }
  
  
  ### Evaluate loss of true parameters (up to a constant)
  n0_vec_list_true = lapply(n0_vec_list_true, function(vec) vec-min(vec))
  loss = eval_loss(edge_time_mat_list = edge_time_mat_list, n0_vec_list = n0_vec_list_true, 
                   clusters_list = clusters_list,  t_vec = t_vec)$loss
  loss_history = c(loss_history, loss)
  

  return(list(n0_vec_list=n0_vec_list_current, gradient_vec_list=gradient_vec_list, 
              n_iter=n_iter, earliest_edge_time_list=earliest_edge_time_list, 
              loss_history=loss_history, accuracy_history=accuracy_history))
  
}







# Test --------------------------------------------------------------------

# res = est_n0_vec_2(edge_time_mat_list = edge_time_mat_list, clusters_list = clusters_list,
#                    t_vec = network$t_vec, step_size = 1e-4, max_iter = 50,
#                    n0_vec_list_true = n0_vec_list_true)
# 
# 
# n0_vec = res$n0_vec_list[[1]]
# gradient_vec = res$gradient_vec_list[[1]]
# gradient_vec * 1e-4
# res$n_iter
# 
# # plot(y=n0_vec*t_unit,x=(network$tau_vec[1:clus_size]-min(network$tau_vec[1:clus_size])))
# plot(y=n0_vec*t_unit,x=(network$tau_vec))
# abline(a=0,b=1,col=2)
# 
# plot(res$loss_history)
# plot(res$accuracy_history)
