t_vec = seq(0,200,1)

empirical_cdf = ecdf(83)(t_vec)
set.seed(100)
event_time_vec = c(runif(30,70,95),rep(Inf,15))
event_time_vec_2 = c(runif(30,70,90),rep(Inf,2))
estimated_cdf = ecdf(event_time_vec)(t_vec)
estimated_cdf_2 = ecdf(event_time_vec_2)(t_vec)

g_cdf = ggplot() +
  geom_line(mapping = aes(x=t_vec, y=empirical_cdf),
            size=1, color='black', alpha=0.6, linetype='solid') +
  geom_line(mapping = aes(x=t_vec, y=estimated_cdf),
            size=1, color='red', alpha=0.6, linetype='dashed') +
  geom_ribbon(aes(x=t_vec, ymin=pmin(empirical_cdf, estimated_cdf),
                  ymax=pmax(empirical_cdf, estimated_cdf)),
              size=0,
              fill="red", alpha=0.5, colour=NA) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 1))
g_cdf

g_cdf_2 = ggplot() +
  geom_line(mapping = aes(x=t_vec, y=empirical_cdf),
            size=1, color='black', alpha=0.6, linetype='solid') +
  geom_line(mapping = aes(x=t_vec, y=estimated_cdf_2),
            size=1, color='blue', alpha=0.6, linetype='dashed') +
  geom_ribbon(aes(x=t_vec, ymin=pmin(empirical_cdf, estimated_cdf_2),
                  ymax=pmax(empirical_cdf, estimated_cdf_2)), 
              fill="blue", alpha=0.5) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 1))
g_cdf_2
  



freq_trun = 6
empirical_pdf = c(hist(83,breaks = t_vec,plot = FALSE)$density,0)
# 
estimated_pdf_fft = fft(c((sum(event_time_vec<Inf)/length(event_time_vec))*hist(event_time_vec,breaks = t_vec,plot = FALSE)$density,0))/length(t_vec)
estimated_pdf_fft_trun = c(head(estimated_pdf_fft, freq_trun+1),
                           rep(0, length(t_vec)-2*freq_trun-1),
                           tail(estimated_pdf_fft, freq_trun))
estimated_pdf = Re(fft(estimated_pdf_fft_trun, inverse = TRUE))
# 
estimated_pdf_fft_2 = fft(c((sum(event_time_vec_2<Inf)/length(event_time_vec_2))*hist(event_time_vec_2,breaks = t_vec,plot = FALSE)$density,0))/length(t_vec)
estimated_pdf_fft_trun_2 = c(head(estimated_pdf_fft_2, freq_trun+1),
                           rep(0, length(t_vec)-2*freq_trun-1),
                           tail(estimated_pdf_fft_2, freq_trun))
estimated_pdf_2 = Re(fft(estimated_pdf_fft_trun_2, inverse = TRUE))
# 
g_pdf = ggplot() +
  geom_line(mapping = aes(x=t_vec, y=empirical_pdf),
            size=1, color='black', alpha=0.6, linetype='solid') +
  geom_line(mapping = aes(x=t_vec, y=estimated_pdf),
            size=1, color='red', alpha=0.6, linetype='solid') +
  geom_line(mapping = aes(x=t_vec, y=estimated_pdf_2),
            size=1, color='blue', alpha=0.6, linetype='solid') +
  coord_cartesian(ylim=c(-0.02,0.08), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 0.5))
g_pdf

sum((empirical_pdf-estimated_pdf)^2)
sum((empirical_pdf-estimated_pdf_2)^2)



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

