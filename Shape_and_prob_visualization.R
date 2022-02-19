library(ggplot2)
t_vec = seq(0,200,1)

empirical_cdf = ecdf(83)(t_vec)
set.seed(10)
event_time_vec = c(runif(30,45,120),rep(Inf,60))
event_time_vec_2 = c(runif(30,45,120),rep(Inf,1))
estimated_cdf = ecdf(event_time_vec)(t_vec)
estimated_cdf_2 = ecdf(event_time_vec_2)(t_vec)
estimated_cdf_3 = pmin(estimated_cdf_2, max(estimated_cdf))
estimated_cdf_4 = c(rep(0,70),head(estimated_cdf_2, length(estimated_cdf_2)-70))
estimated_cdf_5 = pmin(estimated_cdf_4, max(estimated_cdf))

g_cdf = ggplot() +
  geom_line(mapping = aes(x=t_vec, y=estimated_cdf_2),
            size=1, color='black', alpha=0.6, linetype='solid') +
  geom_line(mapping = aes(x=t_vec, y=estimated_cdf),
            size=1, color='black', alpha=0.6, linetype='dashed') +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(estimated_cdf_2, estimated_cdf_3),
                  ymax=pmax(estimated_cdf_2, estimated_cdf_3)),
              size=0,
              fill="red", alpha=0.5, colour=NA) +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(estimated_cdf, estimated_cdf_3),
                  ymax=pmax(estimated_cdf, estimated_cdf_3)),
              size=0,
              fill="blue", alpha=0.5, colour=NA) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 1))
g_cdf


g_cdf_2 = ggplot() +
  geom_line(mapping = aes(x=t_vec, y=estimated_cdf_4),
            size=1, color='black', alpha=0.6, linetype='solid') +
  geom_line(mapping = aes(x=t_vec, y=estimated_cdf),
            size=1, color='black', alpha=0.6, linetype='dashed') +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(estimated_cdf_4, estimated_cdf_5),
                  ymax=pmax(estimated_cdf_4, estimated_cdf_5)),
              size=0,
              fill="red", alpha=0.5, colour=NA) +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(estimated_cdf, estimated_cdf_5),
                  ymax=pmax(estimated_cdf, estimated_cdf_5)),
              size=0,
              fill="blue", alpha=0.5, colour=NA) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 1))
g_cdf_2





pdf(file=paste0("../Results/Plots/Temp/Shape_and_prob/",
                "cdf",
                ".pdf"),
    width = 3, height = 2)
g_cdf
dev.off()

pdf(file=paste0("../Results/Plots/Temp/Shape_and_prob/",
                "cdf_2",
                ".pdf"),
    width = 3, height = 2)
g_cdf_2
dev.off()

pdf(file=paste0("../Results/Plots/Temp/Shape_and_prob/",
                "pdf",
                ".pdf"),
    width = 2, height = 1.3)
g_pdf
dev.off()

