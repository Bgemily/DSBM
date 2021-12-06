
plot_pdf_array_v3 = function(pdf_array_list, 
                             clus_size_vec_1,
                             clus_size_vec_2,
                             t_vec, 
                             y_lim=c(0, 0.04),x_lim=NULL)
{
  
  # library(tidyverse)
  # library(ggplot2)
  
  if (dim(pdf_array_list[[1]])!=dim(pdf_array_list[[2]]))
    stop("Dimension not match: pdf_array_list[[1]] and pdf_array_list[[2]] should have same dimension")
  
  t_unit = t_vec[2] - t_vec[1]
  
  k1 = dim(pdf_array_list[[1]])[1]
  k2 = dim(pdf_array_list[[1]])[2]
  
  big.df = data.frame()
  for (q in 1:k1) {
    for (l in 1:k2) {
      ### Get intensity between cluster (q,l) for L and R spines
      pdf_array_mat = sapply(pdf_array_list, function(pdf_array)pdf_array[q,l,]) # length(t_vec)*ind_subj
      colnames(pdf_array_mat) = paste0('ind_subj', seq(ncol(pdf_array_mat)))
      ### Keep intensity from L if q>l and R if q<l
      if (q>l) {
        pdf_array_mat = pdf_array_mat[,1,drop=FALSE]
      } else if (q<l) {
        pdf_array_mat = pdf_array_mat[,2,drop=FALSE]
      }
      
      pdf_array_mat = cbind(t=t_vec, pdf_array_mat)
      
      tmp.df = pivot_longer(as.data.frame(pdf_array_mat), 
                            cols=starts_with(c('ind_subj')), 
                            names_to = "ind_subj", 
                            values_to = "pdf_val")
      
      tmp.df$clus.pair = paste0("(",q,",",l,")")
      
      tmp.df$fill.col = grDevices::rainbow(300)[floor(tmp.df$t)+1]
      
      big.df = rbind(big.df, tmp.df)
    }
  }
  
  intensity_color_mat = matrix(palette()[1],nrow=k1,ncol=k2)
  diag(intensity_color_mat) = palette()[2:(k1+1)]
  for (q in 1:k1) {
    for (l in q:k2) {
      if(q < l){
        rgb_ave = (1/2)*(col2rgb(intensity_color_mat[q,q])+
                    col2rgb(intensity_color_mat[l,l]) )
        intensity_color_mat[q,l] = rgb(red = rgb_ave[1],
                                       green = rgb_ave[2],
                                       blue = rgb_ave[3], 
                                       maxColorValue = 255)
        intensity_color_mat[l,q] = intensity_color_mat[q,l]
      }
    }
  }
  big.df = big.df %>%
    group_by(ind_subj, clus.pair) %>% 
    mutate(quantile_5 = sum(cumsum(pdf_val)/sum(pdf_val)<(5/100)),
           quantile_95 = sum(cumsum(pdf_val)/sum(pdf_val)<(95/100))) %>%
    ungroup()
  
  g <- big.df %>%
    mutate(ind_subj = as.factor(ind_subj), 
           clus.pair=as.factor(clus.pair)) %>%
    ggplot(aes(x=t, y=pdf_val, group=ind_subj, 
               color=clus.pair,
               # alpha=type_group, 
               # size=type_group, 
               linetype=ind_subj)) +
    geom_line()+
    geom_vline(aes(xintercept = quantile_5, color=clus.pair, linetype=ind_subj),
               alpha = 0.7) +
    geom_vline(aes(xintercept = quantile_95, color=clus.pair, linetype=ind_subj),
               alpha = 0.7) +
    scale_color_manual(values = intensity_color_mat,
                       name=NULL,aesthetics = "color") +
    # scale_alpha_manual(values=c(`1`=0.7, `0`=0.5), name=NULL) +
    # scale_size_manual(values=c(`1`=0.5, `0`=0.5), name=NULL)+
    # scale_linetype_manual(values=c(`1`=1, `0`=2),) + 
    scale_y_continuous(position = "right", limits = y_lim)+
    facet_wrap(~clus.pair) +
    xlab("Time(min)") + ylab(NULL) +
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
    theme(legend.position = 'none') +
    theme(strip.background = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank())  
  
  label_df = big.df %>%
    group_by(clus.pair, ind_subj) %>%
    summarise(conn_prob = round(sum(pdf_val)*(t_vec[2]-t_vec[1]),digits=2))
  clus_size_vec = rep(NA,nrow(label_df))
  clus_size_vec[label_df$clus.pair %in% c("(1,1)","(2,2)","(3,3)","(4,4)") & 
              label_df$ind_subj == 'ind_subj1'] = clus_size_vec_1
  clus_size_vec[label_df$clus.pair %in% c("(1,1)","(2,2)","(3,3)","(4,4)") & 
              label_df$ind_subj == 'ind_subj2'] = clus_size_vec_2
  
  label_df = label_df %>% 
    ungroup() %>%
    mutate(clus_size = clus_size_vec) %>%
    group_by(clus.pair) %>%
    summarise(ind_subj=ind_subj,
              label = paste0('Pr = ', paste0(conn_prob,collapse=","),
                          ifelse(is.na(clus_size),
                                 yes = '', no = paste0('\n', "  = ",paste0(clus_size,collapse=",") )) ) ) %>%
    ungroup()
  
  
  
  g = g + 
    geom_label(data = label_df, x=Inf, y=Inf, 
               mapping = aes(label=label), size=3, color='black',
               vjust = "inward", hjust = "inward",
               show.legend = FALSE, inherit.aes = FALSE)
  
  if(!is.null(x_lim)){
    g = g + xlim(x_lim)
  }
  
  return(list(g=g, big.df=big.df))
  
  # fill curves with colors
  # g + geom_rect(data=big.df, 
  #               mapping=aes(xmin=t,xmax=t+0.3,ymax=value,ymin=0), 
  #               fill=big.df$fill.col,stat="identity")
  
  
}
