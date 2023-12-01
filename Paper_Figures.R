######## Figures for paper #######
library(tidyverse)
library(ggeffects)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(marginaleffects)
library(nlme)

world <- ne_countries(scale = "medium", returnclass = "sf")

## Load year-round data
total.dat <- read.csv("modeling_data.csv")

full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, total.dat)



square <- function(x){
  x^2
} 
# Analyse model
#lstrends(full.model, pairwise ~ urban2, var="abslat")
# significantly negative slope for all, all significantly different from one another
hypothesis_test(full.model, c("abslat", "urban2")) # ggeffects
# all significantly negative, all significantly different from one another (same as above)
#emtrends(full.model, pairwise ~ hemisphere, var="abslat", by="urban2")
# Southern hemisphere is significantly steeper than northern in each urbanization level

# look at means
marginal_means(full.model, variables=c("abslat", "urban2"), cross=TRUE, transform=square)

avg_slopes(full.model, variables="abslat", by="urban2")


hypothesis_test(full.model, terms=c("abslat", "urban2"), test=NULL, scale="response")
hypothesis_test(full.model, terms=c("abslat", "urban2"), scale="response")

## Load seasonal data
dat.season <- read.csv("season_modeling_data.csv")

season.model <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                     precip + log(number_checklists) + elevation, dat.season)
# look at means
season.means <- marginal_means(season.model, variables=c("abslat", "urban2", "season"), cross=TRUE, transform=square)

# Analyse model
#lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Summer")) # all significantly negative, urban and suburban NS
#lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Winter")) # all significantly negative, urban and suburban significantly different
#lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Summer", hemisphere="northern")) # slope overlaps 0 in urb in N hemisphere, suburban and urban NS diff
#lstrends(season.model, pairwise ~ urban2, var=c("abslat"), at=c(season="Summer", hemisphere="southern")) # slope negative for all, suburb and urb NS diff.

# significantly negative slope for all of them
#lstrends(season.model, pairwise ~ urban2, var="abslat", at=c(season="Winter"))

#hypothesis_test(season.model, c("abslat", "urban2", "season", "hemisphere"), test = NULL, scale=exp) # ggeffects

# urban summer in the northern hemisphere the slope overlaps 0
hypothesis_test(season.model, c("abslat", "urban2", "season"), test = NULL)
# also overlaps 0 in full model because of the way it averages (way more points in N hemsiphere)
hypothesis_test(season.model, c("abslat", "urban2", "season"))
# suburban and urban significantly different in the summer but barely (P=0.045), prob won't be robust to thinning



############### Maps of coverage ####################
# Full
full.map<-ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=total.dat, aes(x=long, y=lat, color=urban2, shape=urban2), size=2, alpha=0.5) +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = TRUE) +
  labs(x="Longitude", y="Latitude")+
  # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
  #  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_classic()+
  theme(legend.title=element_blank(), legend.position = c(.7, .3), text=element_text(size=15), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title=element_blank(), axis.line = element_blank())
full.map
ggsave(full.map, file="coverage.map.png", height=12, width=8)

# Seasonal
season.map<-ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=dat.season, aes(x=long, y=lat, color=urban2, shape=urban2), size=1, alpha=0.5) +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = TRUE) +
  labs(x="Longitude", y="Latitude")+
  # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
  #  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_classic()+
  facet_wrap(~season, ncol=1)+
  theme(legend.title=element_blank(),text=element_text(size=15), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title=element_blank(), axis.line = element_blank(), legend.position = "none")
season.map

map <- full.map / season.map
ggsave(season.map, file="coverage.season.map.png", height=8, width=4)






################# LDG PLOTS ##################
square <- function(x){
  x^2
} # make function to square

# plotting with marginal effects
plot_slopes(full.model, variables="abslat", condition=c("urban2", "hemisphere")) # all very much different

# northern hemisphere is steeper

predicted.full<-avg_predictions(full.model, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))


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
#mainLDGplot | marginal.full




ggsave(mainLDGplot, file="LDGMainResults.png", height=5, width=7)
#plot.full.model | thinned.plots # compare the full model with the thinned model

# ok they are not too different except that the one with marginal means is steeper (which makes sense)

## Plot proportional results
proportion.results <- predicted.full %>% select(estimate, urban2, abslat, conf.high, conf.low) %>% 
  pivot_wider(names_from=c("urban2"), values_from = c("estimate", "conf.high", "conf.low")) %>% 
  mutate(proportion.urb = estimate_Urban/estimate_Natural, UCL.urb = conf.high_Urban/conf.low_Natural, 
         LCL.urb = conf.low_Urban/conf.high_Natural,
         proportion.suburb = estimate_Suburban/estimate_Natural, UCL.suburb = conf.high_Suburban/conf.low_Natural, 
         LCL.suburb = conf.low_Suburban/conf.high_Natural)


proportion.plot <- ggplot(proportion.results, aes(x=abslat, y=proportion.urb))+
  geom_line(aes(x=abslat, y=proportion.urb), color="#000000", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.urb, ymin=LCL.urb), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb), color="#CC79A7", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.suburb, ymin=LCL.suburb), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  theme_classic()+
  labs(y="Proportion of natural diversity", x="Absolute Latitude")+
  theme(text=element_text(size=15))

proportion.plot

ggsave(proportion.plot, file="full.proportion.plot.png", height=5, width=6)




######## Plot results of seasonal data
marginal.season <- plot_predictions(season.model, condition=c("abslat", "urban2", "season"), transform=square, points=0.01)

plot_slopes(season.model, variables="abslat", condition=c("urban2", "hemisphere", "season"))

predicted.season <- avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, 
                                         newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                            season = c("Summer", "Winter")))



#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat.season, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
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

#seasonLDGplot / marginal.season
# yes they look the same, great!

ldg.full.results <- mainLDGplot / seasonLDGplot + plot_annotation(tag_levels = "A") + plot_layout(heights = c(2, 1)) +
  labs(tag = "Species Richness") +
  theme(
    plot.tag = element_text(size = 15, angle = 90),
    plot.tag.position = "left"
  )  

ggsave(ldg.full.results2, file="ldg.results.png", height=10, width=8)



## Plot proportion for seasonal data
season.proportion.results <- predicted.season %>% select(estimate, season, urban2, abslat, conf.high, conf.low) %>% 
  pivot_wider(names_from=c("urban2", "season"), values_from = c("estimate", "conf.high", "conf.low")) %>% 
  mutate(proportion.urb.wint = estimate_Urban_Winter/estimate_Natural_Winter, UCL.urb.wint = conf.high_Urban_Winter/conf.low_Natural_Winter, 
         LCL.urb.wint = conf.low_Urban_Winter/conf.high_Natural_Winter,
         proportion.urb.sum = estimate_Urban_Summer/estimate_Natural_Summer, UCL.urb.sum = conf.high_Urban_Summer/conf.low_Natural_Summer, 
         LCL.urb.sum = conf.low_Urban_Summer/conf.high_Natural_Summer,
         proportion.suburb.wint = estimate_Suburban_Winter/estimate_Natural_Winter, UCL.suburb.wint = conf.high_Suburban_Winter/conf.low_Natural_Winter, 
         LCL.suburb.wint = conf.low_Suburban_Winter/conf.high_Natural_Winter, 
         proportion.suburb.sum = estimate_Suburban_Summer/estimate_Natural_Summer, UCL.suburb.sum = conf.high_Suburban_Summer/conf.low_Natural_Summer, 
         LCL.suburb.sum = conf.low_Suburban_Summer/conf.high_Natural_Summer)


proportion.plot.wint <- ggplot(season.proportion.results, aes(x=abslat, y=proportion.urb.wint))+
  geom_line(aes(x=abslat, y=proportion.urb.wint), color="#000000", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.urb.wint, ymin=LCL.urb.wint), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb.wint), color="#CC79A7", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.suburb.wint, ymin=LCL.suburb.wint), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  theme_classic()+
  labs(y="Proportion of natural diversity", x="Absolute Latitude")+
  theme(text=element_text(size=15), axis.title.x=element_blank())

proportion.plot.sum <- ggplot(season.proportion.results, aes(x=abslat, y=proportion.urb.sum))+
  geom_line(aes(x=abslat, y=proportion.urb.sum), color="#000000", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.urb.sum, ymin=LCL.urb.sum), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb.sum), color="#CC79A7", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.suburb.sum, ymin=LCL.suburb.sum), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  theme_classic()+
  labs(y="Proportion of natural diversity", x="Absolute Latitude")+
  theme(text=element_text(size=15))


proportion.plot.season <- proportion.plot.wint / proportion.plot.sum
ggsave(proportion.plot.wint, file="proportion.plot.wint.png", height=5, width=6)
ggsave(proportion.plot.sum, file="proportion.plot.sum.png", height=5, width=6)


### Trying different way of calculating proportion where I subtract the values from the predicted natural value for each point and then plot that 
## Instead of subtracting the model fits 

# Make prediction for every point
# Make prediction for every category at every observed latitude
predicted.Nhemisphere.winter <- predictions(
  season.model,
  type = "response", transform=square,
  # by = "urban2",
  newdata = datagridcf(urban2=c("Natural", "Urban", "Suburban"))) %>% 
  select(estimate, abslat, urban2, hemisphere, season) %>% filter(hemisphere=="northern", season=="Winter") %>% pivot_wider(names_from="urban2", values_from=c("estimate"), values_fn=mean)
# need to do nothern and southern hemisphere seperate because there are overlapping values

predicted.Shemisphere.winter <- predictions(
  season.model,
  type = "response", transform=square,
  # by = "urban2",
  newdata = datagridcf(urban2=c("Natural", "Urban", "Suburban"))) %>% 
  select(estimate, abslat, urban2, hemisphere, season) %>% filter(hemisphere=="southern", season=="Winter") %>% pivot_wider(names_from="urban2", values_from="estimate", values_fn=mean)


predicted.Nhemisphere.summer <- predictions(
  season.model,
  type = "response", transform=square,
  # by = "urban2",
  newdata = datagridcf(urban2=c("Natural", "Urban", "Suburban"))) %>% 
  select(estimate, abslat, urban2, hemisphere, season) %>% filter(hemisphere=="northern", season=="Summer") %>% pivot_wider(names_from="urban2", values_from=c("estimate"), values_fn=mean)
# need to do nothern and southern hemisphere seperate because there are overlapping values

predicted.Shemisphere.summer <- predictions(
  season.model,
  type = "response", transform=square,
  # by = "urban2",
  newdata = datagridcf(urban2=c("Natural", "Urban", "Suburban"))) %>% 
  select(estimate, abslat, urban2, hemisphere, season) %>% filter(hemisphere=="southern", season=="Summer") %>% pivot_wider(names_from="urban2", values_from="estimate", values_fn=mean)



predicted.fulldata <- rbind(predicted.Nhemisphere.winter, predicted.Shemisphere.winter, predicted.Nhemisphere.summer, predicted.Shemisphere.summer)

# merge with observed values
predicted.fulldata.season <- inner_join(dat.season[, c(3, 9, 10, 20, 21)], predicted.fulldata, by=c("abslat", "hemisphere", "season"))

# line for urban sites
resids <- predicted.fulldata.season %>% group_by(urban2) %>% mutate(resids = total_SR/Natural)

# plot
ggplot(resids, aes(x=abslat, y=resids, group=urban2, color=urban2))+
  geom_point(size=0.2, alpha=0.2)+
  facet_grid(~season)+
  geom_smooth(method="lm")













################# Plot thinned data results

## Full
thinned.results <- read.csv("thinned.results.csv")

#thinned.results.summary <- thinned.results %>% group_by(x, group) %>% summarise(mean_x=mean(predicted), max.conf.high = quantile(conf.high, 0.975), min.conf.low = quantile(conf.low, 0.25))
thinned.results.summary <- thinned.results %>% group_by(abslat, urban2) %>% summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))
avg.slope.nat <- thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Natural")
(198-97.8)/7
avg.slope.sub <- thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Suburban")
(148-84)/7
avg.slope.urb <- thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Urban")
(127-82)/7

thinned.plots <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(total.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15))
# this is the plot with the 95% of the confidence intervals
thinned.plots
ggsave(thinned.plots, file="main.thinned.results.png", height=4, width=6)
# compare with full plot
mainLDGplot | thinned.plots
# slopes are the same, just the confidence intervals are different



## Plot proportions
thinned.proportion.results <- thinned.results.summary %>% select(mean_x, urban2, abslat, max.conf.high, min.conf.low) %>% 
  pivot_wider(names_from=c("urban2"), values_from = c("mean_x", "max.conf.high", "min.conf.low")) %>% 
  mutate(proportion.urb = mean_x_Urban/mean_x_Natural, UCL.urb = max.conf.high_Urban/min.conf.low_Natural, LCL.urb = min.conf.low_Urban/max.conf.high_Natural,
           proportion.suburb = mean_x_Suburban/mean_x_Natural, UCL.suburb = max.conf.high_Suburban/min.conf.low_Natural, 
         LCL.suburb = min.conf.low_Suburban/max.conf.high_Natural)
# full proportion results

ggplot(thinned.proportion.results, aes(x=abslat, y=proportion.urb))+
   geom_line(aes(x=abslat, y=proportion.urb), color="#000000", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.urb, ymin=LCL.urb), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb), color="#CC79A7", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.suburb, ymin=LCL.suburb), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  ylim(0,1)+
  theme_classic()+
  labs(y="Proportion of natural diversity", x="Latitude")+
  theme(axis.title.x=element_blank())+
  theme(text=element_text(size=15))







#### Seasonal
seasonal.thinned.results <- read.csv("thinned.seasonal.results.csv")

seasonal.thinned.results.summary <- seasonal.thinned.results %>% group_by(abslat, urban2, season) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))
avg.slope.nat.sum <- seasonal.thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Natural"& season=="Summer")
(129-82.7)/7
avg.slope.nat.wint <- seasonal.thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Natural"& season=="Winter")
(206-15.5)/7
avg.slope.sub.sum <- seasonal.thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Suburban"& season=="Summer")
(81.4-66.3)/7
avg.slope.sub.wint <- seasonal.thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Suburban"& season=="Winter")
(150-21.3)/7
avg.slope.urb.sum <- seasonal.thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Urban" & season=="Summer")
(71.2-62.3)/7
avg.slope.urb.wint <- seasonal.thinned.results.summary %>% filter(abslat%in%c(0,70) & urban2=="Urban"& season=="Winter")
(123-31)/7

# Plot thinned models together
season.thinned.plot <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.season, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(seasonal.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(seasonal.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season, ncol=1)+
  theme(legend.title=element_blank(), legend.position = c(0.2, 0.9), text=element_text(size=15), strip.text=element_blank())
# this is the plot with the 95% of the confidence intervals
season.thinned.plot
ggsave(season.thinned.plot, file="season.thinned.plot.png", height=8, width=6)

ldg.results2 <- thinned.plots / season.thinned.plot + plot_layout(heights = c(2, 1))+ 
  labs(tag = "Species Richness") +
  theme(
    plot.tag = element_text(size = 15, angle = 90),
    plot.tag.position = "left"
  ) 

ldg.results2

ggsave(ldg.results2, file="ldg.thinned.results.png", height=8, width=6)

## Winter thinned plot
dat.wint <- dat.season %>% filter(season=="Winter")
wint.thinned.results.summary <- seasonal.thinned.results.summary %>% filter(season=="Winter")
winter.thinned.plot <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.wint, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(wint.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(wint.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  ylim(0,425)+
 # facet_wrap(~season, ncol=1)+
  theme(legend.title=element_blank(), text=element_text(size=15), strip.text=element_blank(), legend.position="none")

winter.thinned.plot
ggsave(winter.thinned.plot, file="winter.thinned.plot.png", height=4, width=6)


dat.sum <- dat.season %>% filter(season=="Summer")
sum.thinned.results.summary <- seasonal.thinned.results.summary %>% filter(season=="Summer")
summer.thinned.plot <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.sum, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(sum.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(sum.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  ylim(0,425)+
  # facet_wrap(~season, ncol=1)+
  theme(legend.title=element_blank(), text=element_text(size=15), strip.text=element_blank(), legend.position="none")

summer.thinned.plot
ggsave(summer.thinned.plot, file="summer.thinned.plot.png", height=4, width=6)




############################################
### Proportion plots and analyses


############# Full data
### Look at proportion that is lost (expect proportion to be the same)

full.proportion.results <- predicted.full %>% select(estimate, urban2, abslat, conf.high, conf.low) %>% 
  pivot_wider(names_from=c("urban2"), values_from = c("estimate", "conf.high", "conf.low")) %>% 
  mutate(proportion.urb = estimate_Urban/estimate_Natural, UCL.urb = conf.high_Urban/conf.low_Natural, LCL.urb = conf.low_Urban/conf.high_Natural,
         proportion.suburb = estimate_Suburban/estimate_Natural, UCL.suburb = conf.high_Suburban/conf.low_Natural, 
         LCL.suburb = conf.low_Suburban/conf.high_Natural)
# full proportion results

full.proportion.plot <- ggplot(full.proportion.results, aes(x=abslat, y=proportion.urb))+
  geom_line(aes(x=abslat, y=proportion.urb), color="#000000", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.urb, ymin=LCL.urb), alpha=0.2)+
  geom_line(aes(x=abslat, y=proportion.suburb), color="#CC79A7", lwd=1)+
  geom_ribbon(aes(x=abslat, ymax=UCL.suburb, ymin=LCL.suburb), alpha=0.2)+
  geom_hline(yintercept=1, linetype=2, color="#009E73", lwd=1)+
  ylim(0,1)+
  theme_classic()+
  labs(y="Proportion of natural diversity", x="Absolute Latitude")+
  theme(text=element_text(size=15), axis.title.x=element_blank())


ggsave(full.proportion.plot, file="full.proportion.plot.png", height=5, width=6)


# Empty proportion plot (for predictions slides)
proportion.empty <- ggplot(predicted.proportion, aes(x=abslat, y=proportion.urb))+
  # geom_line(aes(x=abslat, y=proportion.urb), color="red")+
  #geom_line(aes(x=abslat, y=proportion.suburb))+
  ylim(0,1)+
  xlim(0,90)+
  theme_classic()+
  labs(y="Proportion of natural diversity", x="Absolute latitude")+
  theme(text=element_text(size=15))
ggsave(proportion.empty, file="proportion.empty.png", height=5, width=6)


### Trying different way of calculating proportion where I subtract the values from the predicted natural value for each point and then plot that 
## Instead of subtracting the model fits 






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










