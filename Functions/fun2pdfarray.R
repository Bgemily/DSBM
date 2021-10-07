
fun2pdfarray = function(true_pdf_fun_list, tau_mat, membership_true, t_vec = seq(0, 50, 0.05)){
  N_clus = length(true_pdf_fun_list)
  pdf_true_array = array(dim = c(N_clus, N_clus, length(t_vec)))
  for (q in 1:N_clus) {
    for (l in 1:N_clus) {
      tau = min(tau_mat[which(membership_true==q), which(membership_true==l)])
      pdf_true_array[q,l,] = sapply(t_vec-tau, true_pdf_fun_list[[q]][[l]])
    }
  }
  return(pdf_true_array)
}