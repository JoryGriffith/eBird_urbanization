######## FRENCH FIGURES #######
library(tidyverse)

##### Empty prediction plot

x <- runif(min=0, max=90, n=150)
y <- runif(min=0, max=300, n=150)

dat <- data.frame(x,y)
# H1

h1.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
#  geom_abline(aes(slope=-2, intercept=300, color="Natural"), linewidth=1.2) +
 # geom_abline(aes(slope=-1.5, intercept=225, color="Suburban"), linewidth=1.2) +
#  geom_abline(aes(slope=-1, intercept=150, color="Urban"), linewidth=1.2) +
  labs(x="Latitude absolue", y="Richesse en espèces") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,400), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
 # scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
  #                   values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank())
h1.plot
ggsave(h1.plot, file="QCBS_figs/blank.plot.french.png", height=4.5, width=6)




########

total.dat <- read.csv("modeling_data.csv")

## Full
thinned.results <- read.csv("thinned.results.csv")

#thinned.results.summary <- thinned.results %>% group_by(x, group) %>% summarise(mean_x=mean(predicted), max.conf.high = quantile(conf.high, 0.975), min.conf.low = quantile(conf.low, 0.25))
thinned.results.summary <- thinned.results %>% group_by(abslat, urban2) %>% summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))

thinned.plots <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(total.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(labels=c("Naturel", "Banlieue", "Urbain"), values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Latitude absolue", y="Richesse en espèces")+
  scale_y_continuous(breaks=c(0,100,200,300,400,500), limits=c(0,500))+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=18))
# this is the plot with the 95% of the confidence intervals
thinned.plots
ggsave(thinned.plots, file="QCBS_figs/main.thinned.results.png", height=6, width=8)
# compare with full plot
mainLDGplot | thinned.plots
# slopes are the same, just the confidence intervals are different
?scale_color_manual




thinned.proportion.results <- thinned.results.summary %>% select(mean_x, urban2, abslat, max.conf.high, min.conf.low) %>% 
  pivot_wider(names_from=c("urban2"), values_from = c("mean_x", "max.conf.high", "min.conf.low")) %>% 
  mutate(proportion.urb = mean_x_Urban/mean_x_Natural, UCL.urb = max.conf.high_Urban/mean_x_Natural, LCL.urb = min.conf.low_Urban/mean_x_Natural,
         proportion.suburb = mean_x_Suburban/mean_x_Natural, UCL.suburb = max.conf.high_Suburban/mean_x_Natural, 
         LCL.suburb = min.conf.low_Suburban/mean_x_Natural)
# full proportion results

thinned.proportion.plot <- ggplot(thinned.proportion.results, aes(x=abslat, y=proportion.urb))+
  geom_line(aes(x=abslat, y=proportion.urb), color="#000000", lwd=1.5)+
  #  geom_ribbon(aes(x=abslat, ymax=UCL.urb, ymin=LCL.urb), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb), color="#CC79A7", lwd=1.5)+
  #  geom_ribbon(aes(x=abslat, ymax=UCL.suburb, ymin=LCL.suburb), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  ylim(0,1)+
  theme_classic()+
  labs(y="Proportion de diversité naturelle", x="Latitude absolue")+
  theme()+
  theme(text=element_text(size=18))

ggsave(thinned.proportion.plot, file="QCBS_figs/thinned.proportion.plot.png", height=5, width=6)




proportion.empty <- ggplot(thinned.proportion.results, aes(x=abslat, y=proportion.urb))+
  # geom_line(aes(x=abslat, y=proportion.urb), color="red")+
  #geom_line(aes(x=abslat, y=proportion.suburb))+
  ylim(0,1)+
  xlim(0,90)+
  theme_classic()+
  labs(y="Proportion de diversité naturelle", x="Latitude absolue")+
  theme(text=element_text(size=18))
ggsave(proportion.empty, file="QCBS_figs/proportion.empty.french.png", height=5, width=6)

dat.season <- read.csv("season_modeling_data.csv")

seasonal.thinned.results <- read.csv("thinned.seasonal.results.csv")
seasonal.thinned.results.summary <- seasonal.thinned.results %>% group_by(abslat, urban2, season) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))
seasonal.thinned.results.summary <- seasonal.thinned.results.summary %>% mutate(season = factor(season, levels=c("Winter", "Summer")))
dat.season <- dat.season %>%  mutate(season = factor(season, levels=c("Winter", "Summer")))
season.thinned.plot <-  
  ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.season, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(seasonal.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(seasonal.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(labels=c("Naturel", "Banlieue", "Urbain"), x="Latitude absolue", y="Richesse en espèces")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(breaks=c(0,100,200,300,400, 500), limits=c(0,500))+
  facet_wrap(~season, ncol=2)+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=14), strip.text=element_blank())
# this is the plot with the 95% of the confidence intervals
season.thinned.plot
ggsave(season.thinned.plot, file="QCBS_figs/season.thinned.plot.png", height=4, width=8)



dat.nat <- dat.season %>% filter(urban2=="Natural")
nat.thinned.results.summary <- seasonal.thinned.results.summary %>% filter(urban2=="Natural")
season.thinned.plot.nat <-  
  ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.nat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(nat.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(nat.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(labels=c("Naturel", "Banlieue", "Urbain"), x="Latitude absolue", y="Richesse en espèces")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(breaks=c(0,100,200,300,400, 500), limits=c(0,500))+
  facet_wrap(~season, ncol=2)+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=14), strip.text=element_blank())
# this is the plot with the 95% of the confidence intervals
season.thinned.plot.nat
ggsave(season.thinned.plot.nat, file="QCBS_figs/season.thinned.plot.nat.png", height=4, width=8)


season.thinned.plot.tall <-  
  ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.season, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(seasonal.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(seasonal.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(labels=c("Naturel", "Banlieue", "Urbain"), x="Latitude absolue", y="Richesse en espèces")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  scale_y_continuous(breaks=c(0,100,200,300,400, 500), limits=c(0,500))+
  facet_wrap(~season, ncol=1)+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=14), strip.text=element_blank())
# this is the plot with the 95% of the confidence intervals
season.thinned.plot.tall
ggsave(season.thinned.plot.tall, file="QCBS_figs/season.thinned.plot.tall.png", height=6, width=4.5)




season.proportion.results <- seasonal.thinned.results.summary %>% select(mean_x, season, urban2, abslat, max.conf.high, min.conf.low) %>% 
  pivot_wider(names_from=c("urban2", "season"), values_from = c("mean_x", "max.conf.high", "min.conf.low")) %>% 
  mutate(proportion.urb.wint = mean_x_Urban_Winter/mean_x_Natural_Winter, UCL.urb.wint = max.conf.high_Urban_Winter/mean_x_Natural_Winter, 
         LCL.urb.wint = min.conf.low_Urban_Winter/mean_x_Natural_Winter,
         proportion.urb.sum = mean_x_Urban_Summer/mean_x_Natural_Summer, UCL.urb.sum = max.conf.high_Urban_Summer/mean_x_Natural_Summer, 
         LCL.urb.sum = min.conf.low_Urban_Summer/mean_x_Natural_Summer,
         proportion.suburb.wint = mean_x_Suburban_Winter/mean_x_Natural_Winter, UCL.suburb.wint = max.conf.high_Suburban_Winter/mean_x_Natural_Winter, 
         LCL.suburb.wint = min.conf.low_Suburban_Winter/mean_x_Natural_Winter, 
         proportion.suburb.sum = mean_x_Suburban_Summer/mean_x_Natural_Summer, UCL.suburb.sum = max.conf.high_Suburban_Summer/mean_x_Natural_Summer, 
         LCL.suburb.sum = min.conf.low_Suburban_Summer/mean_x_Natural_Summer)


proportion.plot.wint <- ggplot(season.proportion.results, aes(x=abslat, y=proportion.urb.wint))+
  geom_line(aes(x=abslat, y=proportion.urb.wint), color="#000000", lwd=1)+
  #  geom_ribbon(aes(x=abslat, ymax=UCL.urb.wint, ymin=LCL.urb.wint), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb.wint), color="#CC79A7", lwd=1)+
  #  geom_ribbon(aes(x=abslat, ymax=UCL.suburb.wint, ymin=LCL.suburb.wint), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  scale_y_continuous(breaks=c(0.5, 0.75, 1, 1.25, 1.5), limits=c(0.5, 1.6))+
  theme_classic()+
  labs(y="Proportion de diversité naturelle", x="Absolute Latitude")+
  theme(text=element_text(size=18), axis.title.x=element_blank())

proportion.plot.sum <- ggplot(season.proportion.results, aes(x=abslat, y=proportion.urb.sum))+
  geom_line(aes(x=abslat, y=proportion.urb.sum), color="#000000", lwd=1)+
  #  geom_ribbon(aes(x=abslat, ymax=UCL.urb.sum, ymin=LCL.urb.sum), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb.sum), color="#CC79A7", lwd=1)+
  # geom_ribbon(aes(x=abslat, ymax=UCL.suburb.sum, ymin=LCL.suburb.sum), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  scale_y_continuous(breaks=c(0.5, 0.75, 1), limits=c(0.5,1))+
  theme_classic()+
  labs(y="Proportion de diversité naturelle", x="Latitude absolue")+
  theme(text=element_text(size=18), axis.title.x=element_blank())


proportion.plot.season <- proportion.plot.wint / proportion.plot.sum
ggsave(proportion.plot.wint, file="QCBS_figs/proportion.plot.wint.png", height=5, width=6)
ggsave(proportion.plot.sum, file="QCBS_figs/proportion.plot.sum.png", height=5, width=6)











