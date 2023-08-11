
### Generate synthetic networks
generate_network = function(N_node, N_clus=3, total_time=200, t_vec=seq(0,total_time,length.out=1000), 
                            clus_size_vec = rep(N_node/N_clus, N_clus),
                            conn_patt_sep=2, conn_prob_mean=0.8, conn_prob_rad=0, 
                                time_shift_mean_vec=rep(10,N_clus), time_shift_rad=min(time_shift_mean_vec), 
                                time_shift_struc=max, 
                                SEED=NULL, const=20)
{
  if(!is.null(SEED)) set.seed(SEED)
  N_subj = 1
  N_node_vec = rep(N_node, N_subj)
  clus_size_mat=matrix(clus_size_vec, byrow = TRUE, nrow=N_subj, ncol=N_clus)
  

# Generate cluster memberships ---------------------------------------
  
  membership_true_list = clus_true_list = vector(mode="list", length=N_subj)
  for (i in 1:N_subj) {
    membership_tmp = rep(1:N_clus, clus_size_mat[i,])
    membership_true_list[[i]] = membership_tmp
    clus_true_list[[i]] = mem2clus(membership_tmp)
  }
  

# Generate time shifts ----------------------------------------------------

  time_shift_list = vector(mode = "list", length = N_subj)
  
  if (length(time_shift_mean_vec) != N_clus) {
    stop("Length of time_shift_mean_vec should be equal to N_clus.")
  }
  
  time_shift_max_vec = time_shift_mean_vec + time_shift_rad
  time_shift_min_vec = time_shift_mean_vec - time_shift_rad
  
  for (m in 1:N_subj) {
    time_shift_tmp = sapply(1:N_clus, function(k)runif(n = clus_size_mat[m,k], 
                                                       min = time_shift_min_vec[k],
                                                       max = time_shift_max_vec[k]) )
    if (is.list(time_shift_tmp)) {
      time_shift_tmp = unlist(time_shift_tmp)
    }
    time_shift_tmp = c(time_shift_tmp)
    
    time_shift_tmp = time_shift_tmp - min(time_shift_tmp)
    
    time_shift_list[[m]] = time_shift_tmp
  }
  


# Generate connecting pattern functions --------------------------------------------

  mean_mat = matrix(nrow = N_clus, ncol = N_clus)
  var_mat = matrix(nrow = N_clus, ncol = N_clus)

  mean_mat[1,2] = mean_mat[2,1] = const*conn_patt_sep
  var_mat[1,2] = var_mat[2,1] = (const^2)/4*conn_patt_sep^(-2)
  
  mean_mat[1,3] = mean_mat[3,1] = const*conn_patt_sep^2
  var_mat[1,3] = var_mat[3,1] = (const^2)/4*conn_patt_sep^(-2)
  
  mean_mat[2,3] = mean_mat[3,2] = const*conn_patt_sep^(1/2) 
  var_mat[2,3] = var_mat[3,2] = (const^2)/4*conn_patt_sep^(-1)
  
  mean_mat[1,1] = const
  var_mat[1,1] = (const^2)/4
  
  mean_mat[2,2] = const*conn_patt_sep^2
  var_mat[2,2] = (const^2)/4*conn_patt_sep
  
  mean_mat[3,3] = const*(conn_patt_sep)^(3/2)
  var_mat[3,3] = (const^2)/4*conn_patt_sep

  conn_prob_min = conn_prob_mean - conn_prob_rad
  conn_prob_max = conn_prob_mean + conn_prob_rad
  if (conn_prob_max > 1 | conn_prob_min < 0) {
    stop("Invalid conn_prob_rad.")
  }
  
  conn_prob_mat = matrix(ncol=N_clus, nrow=N_clus)
  diag(conn_prob_mat) = seq( from=conn_prob_min, to=conn_prob_max, length.out=N_clus )
  conn_prob_mat[1,2] = conn_prob_mat[2,1] = conn_prob_mat[1,1]
  conn_prob_mat[1,3] = conn_prob_mat[3,1] = conn_prob_mat[2,2]
  conn_prob_mat[2,3] = conn_prob_mat[3,2] = conn_prob_mat[3,3]
  
  if (!isSymmetric(conn_prob_mat))
    stop("Constructed conn_prob_mat is not symmetric.")
  

  ### Evaluate true connecting patterns
  pdf_true_array = cdf_true_array = array(dim = c(N_clus, N_clus, length(t_vec)))
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      pdf_true_array[q,l,] = dgamma(x=t_vec, shape=mean_mat[q,l]^2/var_mat[q,l], 
                                    rate=mean_mat[q,l]/var_mat[q,l]) * conn_prob_mat[q,l]
      
      cdf_true_array[q,l,] = pgamma(q=t_vec, shape=mean_mat[q,l]^2/var_mat[q,l], 
                                    rate=mean_mat[q,l]/var_mat[q,l]) * conn_prob_mat[q,l]
    }
  }
  
  

# Generate edge time matrices ---------------------------------------------

  edge_time_mat_list = list()
  for (m in 1:N_subj) {
    N_node_tmp = N_node_vec[[m]]
    
    edge_time_mat_tmp = matrix(Inf, nrow = N_node_tmp, ncol = N_node_tmp)
    for (q in 1:N_clus) {
      for (l in 1:N_clus) {
        samples = ifelse(test = runif(clus_size_mat[m,q]*clus_size_mat[m,l])<=conn_prob_mat[q,l], 
                         yes = rgamma(n = clus_size_mat[m,q]*clus_size_mat[m,l], 
                                      shape=mean_mat[q,l]^2/var_mat[q,l], 
                                      rate=mean_mat[q,l]/var_mat[q,l]), 
                         no = Inf)
        samples = matrix(samples, clus_size_mat[m,q], clus_size_mat[m,l]) 
        edge_time_mat_tmp[clus_true_list[[m]][[q]], clus_true_list[[m]][[l]]] = samples
      }
    }
    edge_time_mat_tmp[lower.tri(edge_time_mat_tmp)] = t(edge_time_mat_tmp)[lower.tri(edge_time_mat_tmp)] # make it symmetric
    
    ### Compute time shift matrix
    time_shift_vec_tmp = time_shift_list[[m]]
    time_shift_mat_tmp = matrix(nrow = N_node_tmp, ncol = N_node_tmp)
    if (is.function(time_shift_struc)){
      for (i in 1:N_node_tmp) {
        for (j in 1:N_node_tmp) {
          time_shift_mat_tmp[i,j] = time_shift_struc(time_shift_vec_tmp[i], time_shift_vec_tmp[j])
        }
      }
    }
    else
      stop("Invalid time_shift_struc. Supposed to be a function.")
    
    
    edge_time_mat_tmp = edge_time_mat_tmp + time_shift_mat_tmp # add time shifts 
    edge_time_mat_tmp[ edge_time_mat_tmp>total_time ] = Inf
    diag(edge_time_mat_tmp) = Inf
    
    edge_time_mat_list[[m]] = edge_time_mat_tmp
  }
  
  


# Output ------------------------------------------------------------------

  
  return(list(edge_time_mat_list=edge_time_mat_list, 
              membership_true_list=membership_true_list, clus_true_list=clus_true_list,
              time_shift_list=time_shift_list, 
              cdf_true_array=cdf_true_array, pdf_true_array=pdf_true_array
              ))
  
}

