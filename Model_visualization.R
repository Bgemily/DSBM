library(ggplot2)

t_vec = seq(0,200,1)
center_pdf_11 = dgamma(x=t_vec, shape=20^2/80, 
                       rate=20/80)
center_pdf_12 = dgamma(x=t_vec, shape=(20*2)^2/(100/2^2), 
                       rate=(20*2)/(100/2^2))
center_pdf_22 = dgamma(x=t_vec, shape=(20*3)^2/(100*2), 
                       rate=(20*3)/(100*2))

pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "f_11",
                ".pdf"),
    width = 2, height = 1.3)
g_11 = ggplot(mapping = aes(x=t_vec, y=center_pdf_11)) +
  geom_line(size=1, color='blue', alpha=0.6) +
  coord_cartesian(ylim=c(0,0.09), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 0.5))
g_11
dev.off()


pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "f_12",
                ".pdf"),
    width = 2, height = 1.3)
g_12 = ggplot(mapping = aes(x=t_vec, y=center_pdf_12)) +
  geom_line(size=1, alpha=0.6) +
  coord_cartesian(ylim=c(0,0.09), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 0.5))
g_12
dev.off()


pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "f_22",
                ".pdf"),
    width = 2, height = 1.3)
g_22 = ggplot(mapping = aes(x=t_vec, y=center_pdf_22)) +
  geom_line(size=1, alpha=0.6) +
  coord_cartesian(ylim=c(0,0.09), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 0.5))
g_22
dev.off()


pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "f_11_annotate",
                ".pdf"),
    width = 2, height = 1.3)
g_11_annotate = g_11 +
  coord_cartesian(ylim=c(-0.02,0.07), xlim=c(0,200)) +
  geom_segment(mapping = aes(x=120,xend=120,y=0,yend=0.004), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=120,xend=120,y=0+0.008,yend=0.004+0.008), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=120,xend=120,y=0+0.016,yend=0.004+0.016), 
               color='orange', size=1, linetype=1, alpha=1 )
g_11_annotate
dev.off()


center_pdf_11_shifted = c(rep(0,120),head(center_pdf_11, length(center_pdf_11)-120))
g_11_shifted = ggplot(mapping = aes(x=t_vec, y=center_pdf_11_shifted)) +
  geom_line(size=1, color=rgb(.5, .25, .5), alpha=0.8) +
  coord_cartesian(ylim=c(-0.01,0.09), xlim=c(0,200)) +
  theme_void()  +
  theme(panel.background=element_rect(colour="gray",size = 0.5))

pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "lambda_ij",
                ".pdf"),
    width = 2, height = 1.3)
g_11_shifted_annotate = g_11_shifted +
  coord_cartesian(ylim=c(-0.02,0.07), xlim=c(0,200)) +
  geom_segment(mapping = aes(x=120,xend=120,y=0,yend=0.004), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=120,xend=120,y=0+0.008,yend=0.004+0.008), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=120,xend=120,y=0+0.016,yend=0.004+0.016), 
               color='orange', size=1, linetype=1, alpha=1 )
g_11_shifted_annotate
dev.off()


pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "f_11_lambda_ij",
                ".pdf"),
    width = 2, height = 1.5)
g_11_shifted_ornot_annotate = g_11 +
  geom_line(mapping = aes(x=t_vec, y=center_pdf_11_shifted),
            size=1, color=rgb(.5, .25, .5), alpha=0.8) +
  coord_cartesian(ylim=c(-0.01,0.09), xlim=c(0,200)) +
  # 
  geom_segment(mapping = aes(x=40,xend=40,y=0,yend=0.004), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=40,xend=40,y=0+0.008,yend=0.004+0.008), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=40,xend=40,y=0+0.016,yend=0.004+0.016), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  # 
  geom_segment(mapping = aes(x=120,xend=120,y=0,yend=0.004), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=120,xend=120,y=0+0.008,yend=0.004+0.008), 
               color='orange', size=1, linetype=1, alpha=1 ) +
  geom_segment(mapping = aes(x=120,xend=120,y=0+0.016,yend=0.004+0.016), 
               color='orange', size=1, linetype=1, alpha=1 )
g_11_shifted_ornot_annotate
dev.off()



pdf(file=paste0("../Results/Plots/Temp/Model_visualization/",
                "dN_ij",
                ".pdf"),
    width = 2, height = 1.5)
g_11_shifted = ggplot(mapping = aes(x=t_vec, y=center_pdf_11_shifted)) +
  geom_line(size=1, color=rgb(.5, .25, .5), alpha=0.3) +
  coord_cartesian(ylim=c(-0.01,0.09), xlim=c(0,200)) +
  theme_void()  +
  theme(panel.background=element_rect(colour="gray",size = 0.5))
g_11_shifted_tij = g_11_shifted +
  coord_cartesian(ylim=c(-0.02,0.09), xlim=c(0,200)) +
  geom_segment(mapping = aes(x=120+15,xend=120+15,y=0,yend=0.058), 
               color=rgb(.5, .25, .5), size=1) +
  geom_segment(mapping = aes(x=0,xend=200,y=0,yend=0),
               color=rgb(.5, .25, .5), size=0.3, 
               arrow = arrow(length = unit(0.05, "npc")) ) +
  theme(panel.background=element_rect(colour="gray",size = 0.5))
g_11_shifted_tij
dev.off()


