
spectral_clustering = function (W, k) 
{
  n = ncol(W)
  S = rowSums(W)
  D = diag(S)
  # L = D - W
  D_sq_inv = diag(S^(-1/2))
  # L = (diag(1, n) - D_sq_inv%*%W%*%D_sq_inv)
  L = (diag(1, n) - diag(S^(-1))%*%W)
  U = (eigen(L)$vectors)[, ((n - k + 1):n)]
  C = cluster::pam(x = U, k = k)
  return(C$clustering)
}