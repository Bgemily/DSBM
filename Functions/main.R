

main = function(case,SEED=NULL, N_clus, N_overclus=N_clus, MaxIter=10, N_subj=1, 
                bw=1, tau_max=5, conn_prob=1, beta=2, alpha=1, const=40,
                tau_struc=max, standardize=FALSE, total_time=50, step_size=0.02, 
                clus_size_vec = c(30,30,30), conv_thres=1e-3)
{
  total_time_rescaled = 50
  total_time = total_time
  ######## V1
  # t_vec = seq(0, total_time_rescaled, 0.05) 
  ####### V2
  t_vec = seq(0, total_time, length.out = 1000)
  #######################
  N_trial = 10
  
  # set.seed(SEED)
  
  if (case==1) network = generate_network1(SEED=SEED,total_time=total_time, tau_max=tau_max, conn_prob=conn_prob, tau_struc=tau_struc)
  if (case==2) network = generate_network2(SEED=SEED,total_time=total_time, 
                                           tau_max=tau_max, conn_prob=conn_prob, beta=beta, alpha=alpha,const = const,
                                           tau_struc=tau_struc, clus_size_vec = clus_size_vec)
  if (case==3) network = generate_network3(SEED=SEED,total_time=total_time)
  if (case==4) network = generate_network4(total_time)
  
  
  network = del_iso_nodes(network)
  edge_time_mat = network$edge_time_mat
  
  # edge_time_mat = edge_time_mat/total_time*total_time_rescaled
  
  res = do_cluster(edge_time_mat = edge_time_mat, N_clus = N_clus, N_overclus = N_overclus, 
                   MaxIter = MaxIter, N_trial = N_trial, t_vec = t_vec, bw = bw, 
                   standardize = standardize, step_size = step_size, conv_thres = conv_thres)
  
  # record ARI
  res$ARI = get_one_ARI(memb_est_vec = clus2mem(res$clusters), 
                    memb_true_vec = network$membership_true)
  # record risk
  true_pdf_array = fun2pdfarray(true_pdf_fun_list = network$true_pdf_fun_list,
                                tau_mat = network$tau_mat*0, 
                                membership_true = network$membership_true,
                                t_vec = network$t_vec)
  t_unit = network$t_vec[2]
  conn_prob = sum(true_pdf_array[1,1,]*t_unit)
  risk = sqrt(sum((true_pdf_array-res$center_pdf_array)^2)*t_unit/N_clus^2)/conn_prob
  res$risk = risk
  
  # return(list(network=network, clus_result=res))
  return(list(clus_result=res[c('ARI','risk','clusters')]))
}




