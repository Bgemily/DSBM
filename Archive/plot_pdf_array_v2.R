
plot_pdf_array_v2 = function(pdf_array_list, pdf_true_array = NULL, 
                             clus_size_vec,
                          t_vec = seq(0, 350, length.out=1000), 
                          y_lim=c(0, 0.04),x_lim=NULL)
{
  
  # library(tidyverse)
  # library(ggplot2)
  
  if (!is.list(pdf_array_list))
    pdf_array_list = list(pdf_array_list)
    
  t_unit = t_vec[2] - t_vec[1]
  
  k1 = dim(pdf_array_list[[1]])[1]
  k2 = dim(pdf_array_list[[1]])[2]
  # par(mfrow = c(k1,k2))
  
  # if(!is.null(pdf_true_array) & k1==k2){ 
  #   if (k1!=dim(pdf_true_array)[1] | k2!=dim(pdf_true_array)[2])
  #     stop("size of pdf_array and pdf_true_array should match.")
  #   
  #   # find permutation for pdf_array
  #   permn_list = lapply(pdf_array_list, function(pdf_array)find_permn(pdf_array, pdf_true_array, t_unit = t_unit, t_vec=t_vec)$permn)
  #   pdf_array_list = mapply(function(pdf_array, permn)pdf_array[permn,permn,], pdf_array_list, permn_list, SIMPLIFY = FALSE)
  # }
  # else if (k1==k2){
  #   permn_list = lapply(pdf_array_list, function(pdf_array)find_permn(pdf_array, pdf_true_array, t_unit = t_unit, t_vec=t_vec)$permn)
  #   pdf_array_list = mapply(function(pdf_array, permn)pdf_array[permn,permn,], pdf_array_list, permn_list, SIMPLIFY = FALSE)
  # }
  # 
  
  big.df = data.frame()
  for (q in 1:k1) {
    for (l in 1:k2) {
      
      pdf_array_mat = sapply(pdf_array_list, function(pdf_array)pdf_array[q,l,]) # length(t_vec)*N_subj
      colnames(pdf_array_mat) = paste0('N_subj', seq(ncol(pdf_array_mat)))
      pdf_array_mat = cbind(t=t_vec, pdf_array_mat, True=c(pdf_true_array[q,l,]))
      
      tmp.df = pivot_longer(as.data.frame(pdf_array_mat), cols=starts_with(c('N_subj','True')), names_to = "N_subj")
      
      tmp.df$clus.pair = paste0("(",q,",",l,")")
      
      tmp.df$fill.col = grDevices::rainbow(300)[floor(tmp.df$t)+1]
      
      big.df = rbind(big.df, tmp.df)
    }
  }
  
  intensity_color_mat = matrix(palette()[1],nrow=k1,ncol=k2)
  diag(intensity_color_mat) = palette()[2:(k1+1)]
  g <- big.df %>%
    mutate(group = as.factor(N_subj), 
           type_group = as.factor((N_subj=='True')*1), 
           clus.pair=as.factor(clus.pair)) %>%
    ggplot(aes(x=t, y=value, group=group, 
               color=clus.pair,
               alpha=type_group, 
               size=type_group, linetype=type_group)) +
    geom_line()+
    # geom_line(color=big.df$fill.col) +
    scale_color_manual(values = intensity_color_mat,
                       name=NULL,aesthetics = "color") +
    scale_alpha_manual(values=c(`1`=0.7, `0`=0.5), name=NULL) +
    scale_size_manual(values=c(`1`=0.5, `0`=0.5), name=NULL)+
    scale_linetype_manual(values=c(`1`=1, `0`=2),) + 
    scale_y_continuous(position = "right", limits = y_lim)+
    facet_wrap(~clus.pair) +
    xlab("Time(min)") + ylab(NULL) +
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
    theme(legend.position = 'none') +
    theme(strip.background = element_blank(),  strip.text.x = element_blank())  
  
  label_df = big.df %>%
    filter(N_subj == "True") %>%
    group_by(clus.pair) %>%
    summarise(conn_prob = round(sum(value)*(t_vec[2]-t_vec[1]),digits=3)) %>%
    mutate(clus_size = ifelse(clus.pair %in% c("(1,1)","(2,2)","(3,3)","(4,4)","(5,5)"), clus_size_vec, NA)) %>%
    mutate(label = paste0('Pr = ',conn_prob,
                          ifelse(is.na(clus_size),
                                 yes = '', no = paste0('\n', "  = ",round(clus_size))) ) )
    
  
  g = g + 
    geom_label(data = label_df, x=Inf, y=Inf, mapping = aes(label=label), size=3, color='black',
              vjust = "inward", hjust = "inward",
              show.legend = FALSE, inherit.aes = FALSE)
  
  if(!is.null(x_lim)){
    g = g + xlim(x_lim)
  }
    
  g
  
  # fill curves with colors
  # g + geom_rect(data=big.df, 
  #               mapping=aes(xmin=t,xmax=t+0.3,ymax=value,ymin=0), 
  #               fill=big.df$fill.col,stat="identity")

  
}
