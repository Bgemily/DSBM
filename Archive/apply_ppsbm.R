

apply_ppsbm = function(case,SEED=NULL, Qmin=3,Qmax=3,
                       tau_max=5, conn_prob=1, beta=2, alpha=1,const=40,
                       tau_struc=NULL, total_time=50, clus_size_vec=c(30,30,30))
{
  # total_time = 50
  
  # set.seed(SEED)
  
  if (case==1) network = generate_network1(SEED=SEED,total_time=total_time, tau_max=tau_max, conn_prob=conn_prob, tau_struc=tau_struc)
  if (case==2) network = generate_network2(SEED=SEED,total_time=total_time, 
                                           tau_max=tau_max, conn_prob=conn_prob, 
                                           beta=beta,alpha=alpha,const=const,
                                           tau_struc=tau_struc, clus_size_vec = clus_size_vec)
  if (case==3) network = generate_network3(SEED=SEED,total_time=total_time)
  if (case==4) network = generate_network4(total_time)
  
  
  network = del_iso_nodes(network)
  edge_time_mat = network$edge_time_mat
  
  library(ppsbm)
  time.seq = numeric(sum(edge_time_mat<Inf))
  type.seq = numeric(sum(edge_time_mat<Inf))
  current_ind = 1
  for (i in 1:nrow(edge_time_mat)) {
    for (j in 1:ncol(edge_time_mat)) {
      if (edge_time_mat[i,j]<Inf){
        time.seq[current_ind] = edge_time_mat[i,j]
        type.seq[current_ind] = convertNodePair(i, j, n = nrow(edge_time_mat), directed = FALSE)
        current_ind = current_ind+1
      }
    }
  }
  data = list(time.seq=time.seq, type.seq=type.seq, Time=total_time)
  Nijk = statistics(data, nrow(edge_time_mat), K=2^6, directed = FALSE)
  
  # res = mainVEM(data=data, n=nrow(edge_time_mat), Qmin=3, directed=FALSE, method="kernel",
  #               d_part=0, n_perturb=0)
  res = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), d_part=5, 
                Qmin=Qmin, Qmax=Qmax, directed=FALSE, method="hist")[[1]]
  res$clusters = mem2clus(apply(res$tau, 2, which.max)) 

  res$ARI = get_one_ARI(memb_est_vec = clus2mem(res$clusters), 
                    memb_true_vec = network$membership_true)

  # return(list(network=network, clus_result=res))
  return(list(clus_result = list(clusters=res$clusters, ARI=res$ARI)))
}

