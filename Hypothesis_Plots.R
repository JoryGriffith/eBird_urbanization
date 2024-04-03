library(terra)
library(sf)
library(tidyverse)
library(spdep)
library(jtools)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggeffects)
library(ggpubr)

x <- runif(min=0, max=90, n=150)
y <- runif(min=0, max=300, n=150)

dat <- data.frame(x,y)
# H1

h1.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-1.5, intercept=300, color="Natural"), linewidth=1.2) +
  geom_abline(aes(slope=-1.5, intercept=225, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-1.5, intercept=150, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(breaks=c(0,25,50,75),expand=c(0,0))+
  scale_y_continuous(breaks=c(0,500), limits=c(0,500), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title=element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank())
h1.plot


h2.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-2, intercept=300, color="Natural"), linewidth=1.2) +
  geom_abline(aes(slope=-1.5, intercept=225, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-1, intercept=150, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
h2.plot

#ggsave(h1.plot, file="blank.plot.png", height=4.5, width=6)
# H2
h3.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-2, intercept=300, linewidth=1.2, color="#009E73") +
  geom_abline(slope=-1.2, intercept=225, linewidth=1.2, color="#CC79A7") +
  geom_abline(slope=-0.5, intercept=150, linewidth=1.2, color="#000000") +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  theme_classic()+
  # annotate("Text", size=5, label="H2: Differing proportions", color="grey20", x=45, y=380)+
  theme(text=element_text(size=15), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
h3.plot


# Winter
winter.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_point(alpha=0)+
  geom_abline(slope=-2.3, intercept=350, linewidth=1.2, color="#009E73") +
  geom_abline(slope=-1.6, intercept=262.5, linewidth=1.2, color="#CC79A7") +
  geom_abline(slope=-.9, intercept=175, linewidth=1.2, color="#000000") +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(breaks=c(0,500), limits=c(0,500), expand=c(0,0))+
  theme_classic()+
  #  annotate("Text", size=5, label="H4: Winter", color="grey20", x=45, y=380)+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title=element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank())
winter.plot

#ggsave(winter.plot, file="committee_meeting_figs/pres_winter_hypothesis_fig.png", height=4, width=6)


# Summer

summer.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-1.3, intercept=275, color="#009E73", linewidth=1.2) +
  geom_abline(slope=-0.8, intercept=206.25, color="#CC79A7", linewidth=1.2) +
  geom_abline(slope=-0.3, intercept=137.5, color="#000000", linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  theme_classic() +
  theme(text=element_text(size=15), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.title=element_blank(), legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank())
# annotate("Text", size=5, label="H4: Summer", color="grey20", x=45, y=380)
summer.plot





# Winter
winter.plot2 <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_point(alpha=0)+
  geom_abline(slope=-3.3, intercept=350, linewidth=1.2, color="#009E73") +
  geom_abline(slope=-2.4, intercept=275, linewidth=1.2, color="#CC79A7") +
  geom_abline(slope=-1.5, intercept=200, linewidth=1.2, color="#000000") +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  theme_classic()+
  #  annotate("Text", size=5, label="H4: Winter", color="grey20", x=45, y=380)+
  theme(text=element_text(size=15), 
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
winter.plot2

#ggsave(winter.plot, file="committee_meeting_figs/pres_winter_hypothesis_fig.png", height=4, width=6)


# Summer

summer.plot2 <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-1.3, intercept=275, color="#009E73", linewidth=1.2) +
  geom_abline(slope=-0.8, intercept=206, color="#CC79A7", linewidth=1.2) +
  geom_abline(slope=-0.3, intercept=137.5, color="#000000", linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  theme_classic() +
  theme(text=element_text(size=15), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.title=element_blank(), legend.position="none", 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
# annotate("Text", size=5, label="H4: Summer", color="grey20", x=45, y=380)
summer.plot2



### Proportion plots
p1.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=-0.002, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-0.002, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())
p1.plot
ggsave(p1.plot, file="hypothesis_figs/proportionH1.plot.png", height=3, width=3)


p2.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=0, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=0, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())
ggsave(p2.plot, file="hypothesis_figs/proportionH2.plot.png", height=3, width=3)


p3.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=0.0027, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=0.003, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())
ggsave(p3.plot, file="hypothesis_figs/proportionH3.plot.png", height=3, width=3)


p4.summer <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=0.0018, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=0.0027, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())
p4.summer
ggsave(p4.summer, file="hypothesis_figs/proportionH4sum.plot.png", height=3, width=3)

p4.winter <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=0.0018, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=0.0027, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())
p4.winter
ggsave(p4.summer, file="hypothesis_figs/proportionH4wint.plot.png", height=3, width=3)

p4.summer <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=0.002, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=0.0025, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())
ggsave(p3.summer, file="hypothesis_figs/proportionH3sum.plot.png", height=3, width=3)

p4.winter <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=0, intercept=1, color="Natural"), linewidth=1.5, lty=2) +
  geom_abline(aes(slope=0.003, intercept=0.75, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=0.0055, intercept=0.5, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Proportion of natural") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0.25,1), expand=c(0,0))+
  theme_classic() +
  # annotate("Text", size=5, label="H1: Same proportion", color="grey20", x=45, y=380)+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  theme(text=element_text(size=15), legend.position="none", legend.title=element_blank(), 
        axis.title.x=element_blank())

ggsave(p4.winter, file="hypothesis_figs/proportionH3wint.plot.png", height=3, width=3)


library(patchwork)
hypothesis.plot <- (h1.plot + h2.plot + h3.plot) / (winter.plot + summer.plot | winter.plot2 + summer.plot2)

#hypothesis.plot <- ggarrange(h1.plot, h2.plot, h3.plot, winter.plot, summer.plot)

require(grid) 
hypothesis.plot2 <- annotate_figure(hypothesis.plot, left = textGrob("Bird Species Richness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                    bottom = textGrob("Absolute Latitude", gp = gpar(cex = 1.3)))
hypothesis.plot2
ggsave(hypothesis.plot, file="hypothesis_figs/hypothesis.plot.png", height=8, width=11)
ggsave(hypothesis.plot2, file="hypothesis_figs/hypothesis.plot.labels.png", height=8, width=11)




