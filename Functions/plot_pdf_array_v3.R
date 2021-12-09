
plot_pdf_array_v3 = function(pdf_array_list, 
                             clus_size_vec_1,
                             clus_size_vec_2,
                             t_vec, 
                             y_lim=c(0, 0.04),xlim=c(0,200))
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
    mutate(pdf_val_quan_5 = pdf_val[t==quantile_5],
           pdf_val_quan_95 = pdf_val[t==quantile_95]) %>%
    ungroup()
  
  ### Get labels
  label_df = big.df %>%
    group_by(clus.pair, ind_subj) %>%
    summarise(conn_prob = round(sum(pdf_val)*(t_vec[2]-t_vec[1]),digits=2))
  
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
    # group_by(clus.pair) %>%
    mutate(ind_subj=ind_subj,
              label = paste0(conn_prob,
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
                   # size=type_group, 
                   linetype=ind_subj)) +
        geom_line(alpha=1)+
        scale_color_manual(values = intensity_color_mat[q,l],
                           name=NULL,aesthetics = "color") +
        scale_linetype_discrete(name='',labels=c('ind_subj1'=label_df[label_df$clus.pair==paste0('(',q,',',l,')')
                                                                      &label_df$ind_subj=='ind_subj1', "label"][[1]],
                                                 'ind_subj2'=label_df[label_df$clus.pair==paste0('(',q,',',l,')')
                                                                      &label_df$ind_subj=='ind_subj2', "label"][[1]]))+
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
              legend.justification = c("right", "top")) +
        guides(color=FALSE, linetype=guide_legend(direction = 'horizontal')) +
        theme(strip.background = element_blank(),  
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.x = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))  + 
        coord_cartesian(xlim = xlim, ylim = y_lim)
      g = g +
        geom_segment(data = big.df %>% filter(clus.pair==paste0('(',q,',',l,')')) %>%
                       group_by(clus.pair, ind_subj) %>% 
                       summarise(quantile_5=unique(quantile_5),
                                 quantile_95=unique(quantile_95)) %>% 
                       ungroup() %>% mutate(clus.pair=as_factor(clus.pair), 
                                            ind_subj=as.factor(ind_subj)),
                     aes(x=quantile_5,xend=quantile_95, 
                         y=y_lim[1], yend=y_lim[1],
                         color=clus.pair, linetype=ind_subj), 
                     inherit.aes = FALSE,
                     arrow=arrow(angle = 0,
                                 length = unit(0.08, "inches"), 
                                 ends='both'),
                     alpha=0.4, size=1) 
      # g = g + 
      #   geom_text(data = label_df[label_df$clus.pair==paste0('(',q,',',l,')'),], 
      #              x=c(150,200), y=y_lim[2],
      #              mapping = aes(label=label),
      #             size=3, color='black',
      #             # position = position_dodge(width = .9),
      #              vjust = "inward", hjust = "inward",
      #             check_overlap = TRUE,
      #             # label.padding = unit(0.15, "lines"),
      #              show.legend = FALSE, inherit.aes = FALSE)
      
      g_list = c(g_list, list(g))
    }
  }
  layout_mat = matrix(NA,3,3)
  layout_mat[upper.tri(layout_mat, diag = TRUE)] = c(1,2,4,3,5,6)
  g = arrangeGrob(grobs=g_list, layout_matrix=layout_mat)
  
  return(list(g=g, big.df=big.df, g_list=g_list))
  
  # fill curves with colors
  # g + geom_rect(data=big.df, 
  #               mapping=aes(xmin=t,xmax=t+0.3,ymax=value,ymin=0), 
  #               fill=big.df$fill.col,stat="identity")
  
  
}
