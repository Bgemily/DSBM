
plot_pdf_array = function(pdf_array_list, pdf_true_array = NULL, 
                          t_vec = seq(0, 350, length.out=1000), 
                          y_lim=c(0, 0.04))
{
  if (!is.list(pdf_array_list))
    pdf_array_list = list(pdf_array_list)
    
  t_unit = t_vec[2] - t_vec[1]
  
  k1 = dim(pdf_array_list[[1]])[1]
  k2 = dim(pdf_array_list[[1]])[2]
  
  big.df = data.frame()
  for (q in 1:k1) {
    for (l in 1:k2) {
      
      pdf_array_mat = sapply(pdf_array_list, function(pdf_array)pdf_array[q,l,]) # length(t_vec)*N_subj
      tmp.df = data.frame(cbind(t=t_vec, true=pdf_true_array[q,l,], pdf_array_mat, mean=rowMeans(pdf_array_mat))) 
      tmp.df = reshape2::melt(tmp.df, id.vars = 't')
      tmp.df$col.group = sapply(as.vector(tmp.df$variable), switch, true="True", mean="Mean", "Estimate")

      tmp.df$clus.pair = paste0("(",q,",",l,")")
      
      tmp.df$fill.col = grDevices::rainbow(300)[floor(tmp.df$t)+1]
      
      big.df = rbind(big.df, tmp.df)
    }
  }
  
  
  library(ggplot2)
  g <- ggplot(big.df, aes(t, value, group=col.group, color=col.group, alpha=col.group, size=col.group, linetype=col.group)) +
    geom_line()+
    scale_color_manual(values = c("Estimate"="black","True"="red", "Mean"="black"), name=NULL,aesthetics = "color") + 
    scale_alpha_manual(values=c("Estimate"=0.5,"True"=0.5, "Mean"=0.7), name=NULL) +
    scale_size_manual(values=c("Estimate"=0.3,"True"=0.3, "Mean"=0.4), name=NULL)+
    scale_linetype_manual(values = c("Estimate"=1,"True"=1, "Mean"=1)) + 
    ylim(y_lim) +
    facet_wrap(~clus.pair) +
    xlab("Time") + ylab(NULL) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(strip.background = element_blank(),  strip.text.x = element_blank())  

  g
  

  
}
