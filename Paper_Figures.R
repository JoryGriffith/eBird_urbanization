######## Figures for paper #######
library(tidyverse)
library(ggeffects)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(marginaleffects)

world <- ne_countries(scale = "medium", returnclass = "sf")

## Load year-round data
total.dat <- read.csv("modeling_data.csv")

full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   BIOME + log(number_checklists) + elevation, total.dat)

# Analyse model
lstrends(full.model, pairwise ~ urban2, var="abslat")
# significantly negative slope for all, all significantly different from one another
hypothesis_test(full.model, c("abslat", "urban2")) # ggeffects
# all significantly negative, all significantly different from one another (same as above)
emtrends(full.model, pairwise ~ hemisphere, var="abslat", by="urban2")
# Southern hemisphere is significantly steeper than northern in each urbanization level


## Load seasonal data
dat.season <- read.csv("season_modeling_data.csv")

season.model <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                     BIOME + log(number_checklists) + elevation, dat.season)

# Analyse model
lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Summer")) # all significantly negative, urban and suburban NS
lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Winter")) # all significantly negative, urban and suburban significantly different
lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Summer", hemisphere="northern")) # slope overlaps 0 in urb in N hemisphere, suburban and urban NS diff
lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Summer", hemisphere="southern")) # slope negative for all, suburb and urb NS diff.

# significantly negative slope for all of them
lstrends(season.model, pairwise ~ urban2, var="abslat", at=c(season="Winter"))

hypothesis_test(season.model, c("abslat", "urban2", "season", "hemisphere"), test = NULL) # ggeffects
# urban summer in the northern hemisphere the slope overlaps 0
hypothesis_test(season.model, c("abslat", "urban2", "season"), test = NULL)
# also overlaps 0 in full model because of the way it averages (way more points in N hemsiphere)
hypothesis_test(season.model, c("abslat", "urban2", "season"))
# suburban and urban significantly different in the summer but barely (P=0.045), prob won't be robust to thinning




############### Maps of coverage ####################
# Full
full.map<-ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=total.dat, aes(x=long, y=lat, color=urban2), size=0.05, alpha=0.3) +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = FALSE) +
  labs(x="Longitude", y="Latitude")+
  # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
  #  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_classic()+
  theme(legend.title=element_blank(), legend.position = c(.75, .3), text=element_text(size=15), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title=element_blank(), axis.line = element_blank())
full.map

# Seasonal
season.map<-ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=dat.season, aes(x=long, y=lat, color=urban2), size=0.05, alpha=0.3) +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = FALSE) +
  labs(x="Longitude", y="Latitude")+
  # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
  #  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title=element_blank(), axis.line = element_blank(), legend.position="none")+
  facet_wrap(~season, ncol=2)
season.map

full.map / season.map






################# LDG PLOTS ##################
square <- function(x){
  x^2
} # make function to square

# plotting with marginal effects
plot_slopes(full.model, variables="abslat", condition=c("urban2", "hemisphere")) # all very much different
# northern hemisphere is steeper

predicted.full<-avg_predictions(full.model, by=c("abslat", "urban2"), transform=square, newdata="mean")
# want to make this over a more limited set of predictors and somehow make not wiggly
#marginal.full<-plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) # need to figure out how to back transform this

mainLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(total.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mainLDGplot
# yay this plot is the same!
mainLDGplot | marginal.full




ggsave(plot.full.model, file="LDGMainResults.png", height=5, width=7)
plot.full.model | thinned.plots # compare the full model with the thinned model

# ok they are not too different except that the one with marginal means is steeper (which makes sense)



######## Plot results of seasonal data
marginal.season <- plot_predictions(season.model, condition=c("abslat", "urban2", "season"), transform=square, points=0.01)

plot_slopes(season.model, variables="abslat", condition=c("urban2", "hemisphere", "season"))

predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(total.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot
# beautiful!

seasonLDGplot / marginal.season
# yes they look the same, great!

ldg.full.results <- plot.full.model / seasonal.results.plot + plot_annotation(tag_levels = "A")

ldg.full.results2 <- wrap_elements(ldg.full.results) +
  labs(tag = "Species Richness") +
  theme(
    plot.tag = element_text(size = 15, angle = 90),
    plot.tag.position = "left"
  )  

ggsave(ldg.full.results2, file="ldg.results.png", height=6, width=8)









################# Plot thinned data results

## Full
thinned.results <- read.csv("thinned.results.csv")

#thinned.results.summary <- thinned.results %>% group_by(x, group) %>% summarise(mean_x=mean(predicted), max.conf.high = quantile(conf.high, 0.975), min.conf.low = quantile(conf.low, 0.25))
thinned.results.summary <- thinned.results %>% group_by(x, group) %>% summarise(mean_x=mean(predicted), max.conf.high = max(conf.high), min.conf.low = min(conf.low))


thinned.plots <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(total.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary, mapping=aes(x=x, y=mean_x, color=group), lwd=1.5)+
  geom_ribbon(thinned.results.summary, mapping=aes(x=x, ymax=max.conf.high, ymin=min.conf.low, group=group), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
# this is the plot with the 95% of the confidence intervals
thinned.plots


#### Seasonal
seasonal.thinned.results <- read.csv("thinned.seasonal.results.csv")

#seasonal.thinned.results.summary <- seasonal.thinned.results %>% group_by(x, group, facet) %>% 
#  summarise(mean_x=mean(predicted), max.conf.high = quantile(conf.high, 0.975), min.conf.low = quantile(conf.low, 0.25))


seasonal.results.plot / season.thinned.plot


# Plot thinned models together
ldg.results <- thinned.plots / season.thinned.plot + plot_annotation(tag_levels = "A")

seasonal.thinned.results.summary <- seasonal.thinned.results %>% group_by(x, group, facet) %>% 
  summarise(mean_x=mean(predicted), max.conf.high = max(conf.high), min.conf.low = min(conf.low))
season.thinned.plot <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(seasonal.thinned.results.summary, mapping=aes(x=x, y=mean_x, color=group), lwd=1.5)+
  geom_ribbon(seasonal.thinned.results.summary, mapping=aes(x=x, ymax=max.conf.high, ymin=min.conf.low, group=group), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~facet)+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=15), axis.title.y=element_blank())
# this is the plot with the 95% of the confidence intervals
season.thinned.plot
ldg.results2 <- wrap_elements(ldg.results) +
  labs(tag = "Species Richness") +
  theme(
    plot.tag = element_text(size = 15, angle = 90),
    plot.tag.position = "left"
  )
ggsave(ldg.results2, file="ldg.thinned.results.png", height=6, width=8)





########## Old plots with emmeans ##########

## Full model

# Plotting with emmeans
#predicted.ggeffect <- ggemmeans(full.model, terms=c("abslat", "urban2"))
#plot.full.model <-
#  plot(predicted.ggeffect, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
#       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
#  theme_classic()+
#  scale_x_continuous(expand=c(0, 0))+
#  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
#  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'), legend.title=element_blank(), legend.position = c(.8, .85), axis.title = element_blank())
#plot.full.model



## Seasonal model
#predicted.ggeffect.season <- ggemmeans(season.model, terms=c("abslat", "urban2", "season")) 
#seasonal.results.plot<-
#  plot(predicted.ggeffect.season, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
#       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
#  theme_classic()+
#  scale_x_continuous(expand=c(0, 0))+
#  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
#  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'), legend.title=element_blank(), legend.position = "none", axis.title.y=element_blank())










