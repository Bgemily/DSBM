
plot_pdf_array = function(pdf_array_list, 
                             clus_size_vec_1,
                             clus_size_vec_2,
                             t_vec, conn_prob_mat_list=NULL,
                          subj1_name="L", subj2_name="R",
                             y_lim=c(0, 0.04),xlim=c(0,200))
{
  
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
      
      pdf_array_mat = cbind(t=t_vec, pdf_array_mat)
      
      tmp.df = pivot_longer(as.data.frame(pdf_array_mat), 
                            cols=starts_with(c('ind_subj')), 
                            names_to = "ind_subj", 
                            values_to = "pdf_val")
      
      tmp.df$clus.pair = paste0("(",q,",",l,")")
      
      tmp.df$fill.col = grDevices::rainbow(300)[floor(tmp.df$t)+1]
      
      if (!is.null(conn_prob_mat_list)) {
        tmp.df$conn_prob[tmp.df$ind_subj=="ind_subj1"] = conn_prob_mat_list[[1]][q,l]
        tmp.df$conn_prob[tmp.df$ind_subj=="ind_subj2"] = conn_prob_mat_list[[2]][q,l]
      } 
      
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
  
  if (is.null(conn_prob_mat_list)) {
    big.df = big.df %>%
      group_by(clus.pair, ind_subj) %>%
      mutate(conn_prob = round(sum(pdf_val)*(t_vec[2]-t_vec[1]),digits=2)) %>%
      ungroup()
  }
  ### Get labels
  label_df = big.df %>%
    group_by(clus.pair, ind_subj) %>%
    summarise(conn_prob = mean(conn_prob))
  
  for (q in 1:k1) {
    for (l in 1:k2) {
      label_df[label_df$clus.pair==paste0('(',q,',',l,')') & 
                 label_df$ind_subj == 'ind_subj1',
               'clus_size'] = paste0("(", clus_size_vec_1[q], "x", clus_size_vec_1[l], ")")
      label_df[label_df$clus.pair==paste0('(',q,',',l,')') & 
                 label_df$ind_subj == 'ind_subj2',
               'clus_size'] = paste0("(", clus_size_vec_2[q], "x", clus_size_vec_2[l], ")")
    }
  }
  
  label_df = label_df %>%
    ungroup() %>%
    mutate(ind_subj=ind_subj,
              label = paste0(ifelse(ind_subj=='ind_subj1',subj1_name,subj2_name), '\n',
                             conn_prob,
                             ifelse(is.na(clus_size),
                                    yes = '', 
                                    no = paste0('\n', clus_size )) ) ) %>%
    ungroup()
  
  
  ### Draw plots
  g_list = list()
  for (q in 1:k1) {
    for (l in q:k2) {
      g <- big.df %>%
        mutate(ind_subj = as.factor(ind_subj), 
               clus.pair=as.factor(clus.pair)) %>%
        filter(clus.pair==paste0('(',q,',',l,')')) %>%
        ggplot(aes(x=t, y=pdf_val, group=ind_subj, 
                   color=clus.pair,
                   size=ind_subj,
                   linetype=ind_subj)) +
        geom_line(alpha=1)+
        scale_color_manual(values = intensity_color_mat[q,l],
                           name=NULL,aesthetics = "color") +
        scale_linetype_discrete(name='',labels=c('ind_subj1'=label_df[label_df$clus.pair==paste0('(',q,',',l,')')
                                                                      &label_df$ind_subj=='ind_subj1', "label"][[1]],
                                                 'ind_subj2'=label_df[label_df$clus.pair==paste0('(',q,',',l,')')
                                                                      &label_df$ind_subj=='ind_subj2', "label"][[1]]))+
        scale_size_manual(name=NULL,values=c(0.5,1))+
        scale_y_continuous(position = "right")+
        facet_wrap(~clus.pair) +
        xlab(NULL) + ylab(NULL) +
        theme_bw() +
        theme(axis.ticks.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(), 
              axis.title.x = element_blank()) +
        theme(legend.position = c(.95, .95),
              legend.background = element_blank(),
              legend.key.size = unit(-1, 'cm'), 
              legend.justification = c("right", "top")) +
        guides(color=FALSE, size=FALSE, linetype=guide_legend(direction = 'horizontal', 
                                                  override.aes = list(size = 0))
               ) +
        theme(strip.text.x = element_blank(),
              strip.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))  + 
        coord_cartesian(xlim = xlim, ylim = y_lim)

      g_list = c(g_list, list(g))
    }
  }
  
  g_list_2 = c(g_list, g_list[2], g_list[3], g_list[5])
  layout_mat = matrix(NA,k1,k2)
  layout_mat[lower.tri(layout_mat, diag = TRUE)] = c(1:(k1*(k1+1)/2))
  layout_mat[upper.tri(layout_mat)] = c(7,8,9)
  g = arrangeGrob(grobs=g_list_2, layout_matrix=layout_mat)
  
  return(list(g=g, big.df=big.df, g_list=g_list))
  

}
