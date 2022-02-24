library(ggplot2)
t_vec = seq(0,200,1)

set.seed(10)
event_time_vec_1 = c(runif(30,20,140),rep(Inf,60))
event_time_vec_2 = c(runif(300,45,120),rep(Inf,50))
event_time_vec_3 = c(runif(300,95,110),rep(Inf,100))

Nbar = ecdf(event_time_vec_1)(t_vec)
# F_true = ecdf(event_time_vec_2)(t_vec)
F_true = pnorm(t_vec,80,20)*(6/7)
# arbitrF = ecdf(event_time_vec_3)(t_vec)
arbitrF = pnorm(t_vec, 100,3)*(1/3+0.1)
estimated_cdf_3 = pmin(F_true, max(Nbar))

F_shift_right = c(rep(0,60),head(F_true, length(F_true)-60))
estimated_cdf_5 = pmin(F_shift_right, max(Nbar))

ind_red = which(Nbar==max(Nbar) & F_true==max(F_true))
ind_blue = setdiff(1:length(t_vec), ind_red)
ind_blue = c(1,ind_blue+1)


g_cdf_1 = ggplot() +
  geom_ribbon(aes(x=t_vec,
                  ymin=pmin(Nbar, F_true),
                  ymax=pmax(Nbar, F_true)),
              size=0,
              fill="blue", alpha=0.4, colour=NA) +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(arbitrF, Nbar),
                  ymax=pmax(arbitrF, Nbar)),
              size=0,
              fill="red", alpha=0.4, colour=NA) +
  geom_line(mapping = aes(x=t_vec, y=arbitrF),
            size=1, color='red', alpha=0.8, linetype='dashed') +
  geom_line(mapping = aes(x=t_vec, y=F_true),
            size=1, color='blue', alpha=0.8, linetype='dashed') +
  geom_line(mapping = aes(x=t_vec, y=Nbar),
            size=1, color='black', alpha=0.8, linetype='solid') +
  coord_cartesian(ylim=c(0,1), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 1))
g_cdf_1


g_cdf_2 = ggplot() +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(Nbar, F_true),
                  ymax=pmax(Nbar, F_true)),
              size=0,
              fill="blue", alpha=0.4, colour=NA) +
  geom_ribbon(aes(x=t_vec, 
                  ymin=pmin(Nbar, F_shift_right),
                  ymax=pmax(Nbar, F_shift_right)),
              size=0,
              fill="red", alpha=0.4, colour=NA) +
  geom_line(mapping = aes(x=t_vec, y=F_shift_right),
            size=1, color='red', alpha=0.8, linetype='dashed') +
  geom_line(mapping = aes(x=t_vec, y=F_true),
            size=1, color='blue', alpha=0.8, linetype='dashed') +
  geom_line(mapping = aes(x=t_vec, y=Nbar),
            size=1, color='black', alpha=0.8, linetype='solid') +
  coord_cartesian(ylim=c(0,1), xlim=c(0,200)) +
  theme_void() +
  theme(panel.background=element_rect(colour="gray",size = 1))
g_cdf_2





pdf(file=paste0("../Results/Plots/Temp/Shape_and_prob/",
                "cdf",
                ".pdf"),
    width = 4, height = 2)
g_cdf_1
dev.off()

pdf(file=paste0("../Results/Plots/Temp/Shape_and_prob/",
                "cdf_2",
                ".pdf"),
    width = 4, height = 2)
g_cdf_2
dev.off()


