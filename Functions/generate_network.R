

# case 1 --------------------------------------------------------------------

generate_network1 = function(SEED=NULL,total_time=50, tau_max=c(10,10), conn_prob=1, tau_struc=NULL)
{
  
  if(!is.null(SEED)) set.seed(SEED)

    
  t_vec = seq(0, total_time, 0.05)
  dist_thres = Inf
  N_clus = 2
  clus_size_vec = c(30,30)
  L_mat = matrix(c(1,0,1,1), nrow=N_clus, ncol=N_clus)
  
  if(length(tau_max)!=N_clus)
    tau_max = rep(tau_max,time=N_clus)
  
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)
  
  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_vec = c(runif(clus_size_vec[1], 0, tau_max[1]), runif(clus_size_vec[2], 0, tau_max[2]))
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x,20,2)*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,20,2),Inf)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  # pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=1, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  # pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  # pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1));  
  # pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  # pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dgamma(x, 2,0.02)*conn_prob,
  # random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,2,0.02),Inf))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dgamma(x, 40,2)*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40,2),Inf))
  pdfNrdsamp_fun_list[[1]][[1]] = pdfNrdsamp_fun_list[[2]][[2]]
  
  # pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  if(is.null(tau_struc)){
    for (q in 1:N_clus) {
      for (l in 1:N_clus) {
        clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
        clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
        if (q==l){
          tmp = mapply(max, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
          tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
        }
        else if (L_mat[q,l]==1){
          tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
          tau_mat[clus_l, clus_q] = t(tau_mat[clus_q, clus_l])
        }
      }
    }
    
  }
  else if (is.function(tau_struc)){
    for (i in 1:N_node) {
      for (j in 1:N_node) {
        tau_mat[i,j] = tau_struc(tau_vec[i], tau_vec[j])
      }
    }
  }
  else
    stop("Invalid tau_struc.")

  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
  
}



# case 2 ---------------------------------------


generate_network2 = function(total_time=50, tau_max=c(10,10,10), conn_prob=1, beta=2, alpha=1,
                             tau_struc=max, SEED=NULL, clus_size_vec = c(30,30,30), const=40)
{
  if(!is.null(SEED)) set.seed(SEED)
  
  ######### V1
  # t_vec = seq(0, total_time, 0.05)
  ##########  V2
  t_vec = seq(0, total_time, length.out = 1000)
  ###########################
  dist_thres = 6
  N_clus = 3
  # clus_size_vec = c(30,30,30)
  L_mat = matrix(c(1,0,0,1,1,0,1,1,1), nrow=N_clus, ncol=N_clus)
  
  if(length(tau_max)!=N_clus)
    tau_max = rep(tau_max,time=N_clus)
  
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)
  
  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  ### V1
  tau_vec = c(runif(clus_size_vec[1], 0, tau_max[1]), runif(clus_size_vec[2], 0, tau_max[2]), runif(clus_size_vec[3], 0, tau_max[3]))
  ### V2
  # tau_vec = c(tau_max[1]*rbinom(clus_size_vec[1], 1, 0.5),
  #             tau_max[2]*rbinom(clus_size_vec[2], 1, 0.5),
  #             tau_max[3]*rbinom(clus_size_vec[3], 1, 0.5))
  ##########
  
  
  ################ V1
  # # design true connecting patterns (N_clus*N_clus) and sampling functions
  # pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  # pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, 40*1,1/beta)*conn_prob, 
  #                                      random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40*1,1/beta),Inf)); 
  # pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  # pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, 40/beta,1)*conn_prob, 
  #                                      random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40/beta,1),Inf));
  # pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]]
  # pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dgamma(x, 40*1,1*beta)*conn_prob,
  #                                      random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40*1,1*beta),Inf));
  # pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  # 
  # pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dgamma(x, 40*1, 1*1)*conn_prob,
  #                                      random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40*1, 1*1),Inf))
  # pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dgamma(x, 40, 1*sqrt(beta))*conn_prob,
  #                                      random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40*beta, 1*1),Inf))
  # pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dgamma(x, 40*1, 1/sqrt(beta))*conn_prob, 
  #                                      random = function(n) ifelse(runif(n)<=conn_prob,rgamma(n,40*1, 1/sqrt(beta)),Inf))
  # # pdfNrdsamp_fun_list[[1]][[1]] = pdfNrdsamp_fun_list[[3]][[3]]
  # # pdfNrdsamp_fun_list[[2]][[2]] = pdfNrdsamp_fun_list[[1]][[1]]
  # 
  ############## V2
  mean_mat = matrix(nrow = N_clus, ncol = N_clus)
  var_mat = matrix(nrow = N_clus, ncol = N_clus)
  
    ################ V1
    # mean_mat[1,2] = mean_mat[2,1] = 40*beta
    # var_mat[1,2] = var_mat[2,1] = 40*beta^2
    # 
    # mean_mat[1,3] = mean_mat[3,1] = 40/beta
    # var_mat[1,3] = var_mat[3,1] = 40/beta
    # 
    # mean_mat[2,3] = mean_mat[3,2] = 40/beta
    # var_mat[2,3] = var_mat[3,2] = 40/beta^2
    # 
    # mean_mat[1,1] = 40*1
    # var_mat[1,1] = 40*1
    # 
    # mean_mat[2,2] = 40*(beta)
    # var_mat[2,2] = 40*beta
    # 
    # mean_mat[3,3] = 40*sqrt(beta)
    # var_mat[3,3] = 40*beta
    ############### V2
    mean_mat[1,2] = mean_mat[2,1] = const*beta
    var_mat[1,2] = var_mat[2,1] = const*alpha*beta^2
    
    mean_mat[1,3] = mean_mat[3,1] = const/beta^2
    var_mat[1,3] = var_mat[3,1] = const*alpha/beta
    
    mean_mat[2,3] = mean_mat[3,2] = const/beta 
    var_mat[2,3] = var_mat[3,2] = const*alpha/beta^2 
    
    mean_mat[1,1] = const
    var_mat[1,1] = const*alpha
    
    mean_mat[2,2] = const*beta^2
    var_mat[2,2] = const*alpha*beta
    
    mean_mat[3,3] = const*sqrt(beta)
    var_mat[3,3] = const*alpha*beta
    ####################################
  
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=mean_mat[1,2]^2/var_mat[1,2], 
                                                                rate=mean_mat[1,2]/var_mat[1,2])*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob, rgamma(n, shape=mean_mat[1,2]^2/var_mat[1,2], 
                                                                                              rate=mean_mat[1,2]/var_mat[1,2]), Inf)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=mean_mat[1,3]^2/var_mat[1,3], 
                                                                rate=mean_mat[1,3]/var_mat[1,3])*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob, rgamma(n, shape=mean_mat[1,3]^2/var_mat[1,3], 
                                                                                               rate=mean_mat[1,3]/var_mat[1,3]), Inf));
  pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]]
  pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dgamma(x, shape=mean_mat[2,3]^2/var_mat[2,3], 
                                                                rate=mean_mat[2,3]/var_mat[2,3])*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob, rgamma(n, shape=mean_mat[2,3]^2/var_mat[2,3], 
                                                                                               rate=mean_mat[2,3]/var_mat[2,3]), Inf));
  pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dgamma(x, shape=mean_mat[1,1]^2/var_mat[1,1], 
                                                                rate=mean_mat[1,1]/var_mat[1,1])*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob, rgamma(n, shape=mean_mat[1,1]^2/var_mat[1,1], 
                                                                                               rate=mean_mat[1,1]/var_mat[1,1]), Inf))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dgamma(x, shape=mean_mat[2,2]^2/var_mat[2,2], 
                                                                rate=mean_mat[2,2]/var_mat[2,2])*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob, rgamma(n, shape=mean_mat[2,2]^2/var_mat[2,2], 
                                                                                               rate=mean_mat[2,2]/var_mat[2,2]), Inf))
  pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dgamma(x, shape=mean_mat[3,3]^2/var_mat[3,3], 
                                                                rate=mean_mat[3,3]/var_mat[3,3])*conn_prob, 
                                       random = function(n) ifelse(runif(n)<=conn_prob, rgamma(n, shape=mean_mat[3,3]^2/var_mat[3,3], 
                                                                                               rate=mean_mat[3,3]/var_mat[3,3]), Inf))
  
  
  ######################################
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  if(is.null(tau_struc)){
    print("no input tau_struc.")
    for (q in 1:N_clus) {
      for (l in 1:N_clus) {
        clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
        clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
        if (q==l){
          tmp = mapply(max, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
          tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
        }
        else if (L_mat[q,l]==1){
          tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
          tau_mat[clus_l, clus_q] = t(tau_mat[clus_q, clus_l])
        }
      }
    }
    
  }
  else if (is.function(tau_struc)){
    for (i in 1:N_node) {
      for (j in 1:N_node) {
        tau_mat[i,j] = tau_struc(tau_vec[i], tau_vec[j])
      }
    }
  }
  else
    stop("Invalid tau_struc.")
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list = pdfNrdsamp_fun_list, tau_mat = tau_mat, clus_true = clus_true, 
                                      clus_size_vec = clus_size_vec, pairwise_dist = pairwise_dist, dist_thres = dist_thres,
                                      t_max = total_time)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
  
}





# case 3 ------------------------------------------------------------------


generate_network3 = function(SEED=0, total_time=50)
{
  # set.seed(SEED)
  # 
  t_vec = seq(0, total_time, 0.05)
  dist_thres = 2
  N_clus = 3
  clus_size_vec = c(30,30,30)
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)
  
  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_vec = c(runif(clus_size_vec[1], 25, 30), runif(clus_size_vec[2],0,10), rep(0, clus_size_vec[3]))
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=10, rate=2), random = function(n) rgamma(n,10,2)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=10, rate=1), random = function(n) rgamma(n,10,1)); 
  pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dgamma(x, shape=20, rate=4), random = function(n) rgamma(n,20,4));  
  pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
      clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
      if (q==l){
        tmp = mapply(max, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
        # tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
        tau_mat[clus_q, clus_l] = 0
      }
      else if (q<l)
        tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
      
      else
        tau_mat[clus_q, clus_l] = t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q))
    }
  }
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
  
}

# case 4 --------------------------------------------------------------------


generate_network4 = function(SEED=0, total_time=50)
{
  # set.seed(SEED)
  
  t_vec = seq(0, total_time, 0.05)
  dist_thres = 2
  N_clus = 3
  clus_size_vec = c(50,50,50)
  N_node = sum(clus_size_vec)
  
  membership_true = rep(1:N_clus, clus_size_vec)
  clus_true = mem2clus(membership_true)

  centers_x = runif(N_node, 0, 1)
  centers_y = runif(N_node, 0, 6)
  node_loc_mat = cbind(centers_x, centers_y)
  
  pairwise_dist = rdist::pdist(node_loc_mat)
  
  tau_vec = c(runif(clus_size_vec[1], 40, 42), runif(clus_size_vec[2],0,2), rep(0, clus_size_vec[3]))
  
  
  # design true connecting patterns (N_clus*N_clus) and sampling functions
  pdfNrdsamp_fun_list = apply(matrix(nrow=N_clus, ncol=N_clus), 1, as.list)
  pdfNrdsamp_fun_list[[1]][[2]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  pdfNrdsamp_fun_list[[2]][[1]] = pdfNrdsamp_fun_list[[1]][[2]]
  pdfNrdsamp_fun_list[[1]][[3]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5)); 
  pdfNrdsamp_fun_list[[3]][[1]] = pdfNrdsamp_fun_list[[1]][[3]] 
  pdfNrdsamp_fun_list[[2]][[3]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1));  
  pdfNrdsamp_fun_list[[3]][[2]] = pdfNrdsamp_fun_list[[2]][[3]]
  
  pdfNrdsamp_fun_list[[1]][[1]] = list(pdf = function(x) dgamma(x, shape=2, rate=.5), random = function(n) rgamma(n,2,0.5))
  pdfNrdsamp_fun_list[[2]][[2]] = list(pdf = function(x) dnorm(x, 5, 1), random = function(n) rnorm(n,5,1))
  pdfNrdsamp_fun_list[[3]][[3]] = list(pdf = function(x) dunif(x, 0, 1*total_time), random = function(n) runif(n,0,1*total_time))
  
  
  # extract true pdf functions, N_clus*N_clus
  true_pdf_fun_list = lapply(pdfNrdsamp_fun_list, lapply, '[[', 'pdf') 
  
  
  # obtain tau_mat from tau_vec, decided by who dominates who
  tau_mat = matrix(nrow = N_node, ncol = N_node)
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      clus_q = clus_true[[q]]; clus_l = clus_true[[l]]
      clus_size_q = clus_size_vec[q]; clus_size_l = clus_size_vec[l]
      if (q==l){
        tmp = mapply(min, matrix(tau_vec[clus_q], clus_size_q, clus_size_l), t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q)))
        tau_mat[clus_q, clus_l] = matrix(tmp, nrow = clus_size_q)
      }
      else if (q<l)
        tau_mat[clus_q, clus_l] = matrix(tau_vec[clus_q], clus_size_q, clus_size_l) 
      
      else
        tau_mat[clus_q, clus_l] = t(matrix(tau_vec[clus_l], clus_size_l, clus_size_q))
    }
  }
  
  
  # generate edge_time_mat
  edge_time_mat = gener_edge_time_mat(pdfNrdsamp_fun_list, tau_mat, clus_true, clus_size_vec, pairwise_dist, dist_thres)
  
  
  return(list(edge_time_mat=edge_time_mat, node_loc_mat=node_loc_mat, tau_vec=tau_vec, tau_mat=tau_mat, 
              true_pdf_fun_list=true_pdf_fun_list, membership_true=membership_true,
              t_vec = t_vec, dist_thres=dist_thres, pairwise_dist=pairwise_dist))
}



