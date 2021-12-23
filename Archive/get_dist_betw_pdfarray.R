
# compute distance between two pdf_array
get_dist_betw_pdfarray = function(pdf_array_1, pdf_array_2, symmetric=FALSE, weights=NULL, t_unit=0.05, t_vec = seq(0, 50, 0.05), use_shift_inv_dist=TRUE, pp=TRUE){ 
  
  # step_size = 4e-4 * max(t_vec)
  t_unit = t_vec[2] - t_vec[1]
  step_size = 0.02
  
  
  if(symmetric){
    ### assuming the connecting pattern matrix is symmetric and square. Used for brute-force search of permutation.
    dist = 0
    n0_BetwSubj_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    for (i in 1:(dim(pdf_array_1)[1])) {
      for (j in (i):dim(pdf_array_1)[2]) {
        if(use_shift_inv_dist){
          res = align_pdf_gd(pdf1 = pdf_array_1[i,j,], pdf2 = pdf_array_2[i,j,], t_unit=t_unit, step_size = step_size)
          
          dist = dist + res$dist_min
          n0_BetwSubj_mat[i,j] = res$n0
          n0_BetwSubj_mat[j,i] = n0_BetwSubj_mat[i,j]
        }
        else{
          dist = dist + norm(as.matrix(pdf_array_1[i,j, ]-pdf_array_2[i,j, ]), type="F")
        }
      }
    }
    return(list(dist = dist, n0_mat = n0_BetwSubj_mat))
  }
  
  else{
    dist_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    n0_BetwSubj_mat = matrix(nrow=dim(pdf_array_1)[1], ncol=dim(pdf_array_1)[2])
    for (i in 1:(dim(pdf_array_1)[1])) {
      for (j in (1):dim(pdf_array_1)[2]) {
        if(use_shift_inv_dist){
          if (pp){
            res = align_pdf_gd(pdf1 = pdf_array_1[i,j,], pdf2 = pdf_array_2[i,j,], t_unit=t_unit, step_size = step_size)
            dist_mat[i,j] = res$dist_min
            n0_BetwSubj_mat[i,j] = res$n0
          }
          else{ # cdf
            res = align_curves_gd(f1 = pdf_array_1[i,j,], f2 = pdf_array_2[i,j,], n0 = 0, pp=FALSE, t_unit=t_unit, step_size = step_size)
            dist_mat[i,j] = res$dist_min
            n0_BetwSubj_mat[i,j] = res$n0
          }
          
        }
        else{
          dist_mat[i,j] = norm(as.matrix(pdf_array_1[i,j, ]-pdf_array_2[i,j, ]), type="F")
        }
      }
    }
    ############ V1
    # if(is.null(weights)) dist = mean(dist_mat)
    # else dist = sum(dist_mat*weights)
    ############ V2
    if(is.null(weights)) dist = sqrt(mean(dist_mat^2))
    else dist = sqrt(sum(dist_mat^2*weights))
    ##########################

    return(list(dist = dist, n0_mat = n0_BetwSubj_mat))
  }
  
}

