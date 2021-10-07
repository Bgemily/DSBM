
# Visualize summary statistics (median & quantile)  ---------------------------------------------------------------------


plot_jitter_boxplot = function(data, ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)), 
                               ylab=NULL, xlab=NULL){  
  library(ggplot2)
  library(tidyverse)
  library(viridis)
  ylim = ylim
  data = pivot_longer(data = data, cols = everything())
  group_names_complete = data$name
  group_names_num = gsub(x=group_names_complete, pattern = "[a-z_]",replacement = "")
  group_names_num = as.factor(as.numeric(group_names_num))
  common_string = gsub(group_names_complete, pattern="[0-9.]", replacement="")[1]
  common_string = substr(common_string, start=1, stop = nchar(common_string)-1)
  if (is.null(xlab)) {
    xlab = common_string
  }
  data$name = group_names_num
  sample_size = data %>% group_by(name) %>% summarize(num=n())
  data %>%
    left_join(sample_size) %>%
    # mutate(myaxis = paste0(name)) %>%
    ggplot( aes(x=name, y=value, group=name, fill=name)) +
    geom_violin(width=1, alpha=0.3) +
    geom_boxplot(width=0.3, alpha=0.7) +
    scale_fill_viridis(discrete = TRUE) +
    ylim(ylim)+
    theme_light() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("") +
    ylab(ylab) + 
    xlab(xlab)
  # coord_flip()
  
}


### "data": N_trial*N_sim_setting
plot_pointrange = function(data, ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)), 
                           ylab=NULL, xlab=NULL){  
  library(ggplot2)
  library(tidyverse)
  library(viridis)
  ylim = ylim
  data = pivot_longer(data = data, cols = everything())
  group_names_complete = data$name
  group_names_num_tmp = gsub(x=group_names_complete, pattern = "[A-Za-z_]",replacement = "")
  group_names_num = as.factor(as.numeric(group_names_num_tmp))
  if (is.na(group_names_num[1])) {
    group_names_num = as.factor((group_names_num_tmp))
  }
  common_string = gsub(group_names_complete, pattern="[0-9.]", replacement="")[1]
  common_string = substr(common_string, start=1, stop = nchar(common_string)-1)
  if (is.null(xlab)) {
    xlab = common_string
  }
  data$name = group_names_num

  data %>%
    ggplot( aes(x=name, y=value)) +
    stat_summary(aes(group=1), 
                 fun.min = function(x)quantile(x,0.25),
                 fun.max = function(x)quantile(x,0.75),
                 fun = median) +
    stat_summary(aes(group=1), color="darkgrey",
                 geom="line",
                 fun = "median") +
    coord_cartesian(ylim=ylim)+
    theme_light() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("") +
    ylab(ylab) + 
    xlab(xlab)
  # coord_flip()
  
}


plot_ARI_compr = function(ARI_our, ARI_ppsbm, ylim=c(0,1), continuous=FALSE, reverse=FALSE, n.breaks=10){
  data_tmp=rbind(reshape2::melt(ARI_our), reshape2::melt(ARI_ppsbm))
  method = rep(c("our method","ppsbm"),times=c(nrow(reshape2::melt(ARI_our)), nrow(reshape2::melt(ARI_ppsbm))))
  data_tmp = cbind(data_tmp, method)
  colnames(data_tmp) = c("V","ARI","method")
  
  if(continuous){
    data_tmp$V <- as.numeric(as.character(data_tmp$V))
    width = (max(data_tmp$V)-min(data_tmp$V))/50
  }
  else{
    width = 0.1
  }
  
  # grouped boxplot
  g<- ggplot(data_tmp, aes(x=V, y=ARI, group=method, col=method)) + 
    # geom_boxplot()+
    stat_summary(geom = "line",fun=median)+
    stat_summary(geom = "point",fun=median)+
    stat_summary(geom = "errorbar",width=width,
                 fun.min = function(z) { quantile(z,0.25) },
                 fun.max = function(z) { quantile(z,0.75) }
    )+
    coord_cartesian(ylim = ylim)+
    theme_light() +
    theme(
      axis.title = element_blank(),
      # legend.position="none",
      plot.title = element_text(size=11)
    ) #+
    # scale_x_discrete(guide = guide_axis(n.dodge=2, check.overlap = T))
  
  if(reverse&continuous)
    g + scale_x_reverse(breaks=unique(data_tmp$V))
  else if(continuous)
    g+scale_x_continuous(breaks=unique(data_tmp$V))
  else
    g
}



