################ These are the figures for my committee meeting #########
library(terra)
library(sf)
library(tidyverse)
library(spdep)
library(jtools)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggeffects)
library(spatialreg)
library(ggpubr)

############## Figs for document ################

############## HYPOTHESIS FIGS ###################
# Fake linear regression
x <- runif(min=0, max=90, n=150)
y <- runif(min=0, max=500, n=150)

dat <- data.frame(x,y)
# H1

h1.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-2.2, intercept=400, color="Natural"), linewidth=1.2) +
  geom_abline(aes(slope=-2.2, intercept=320, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-2.2, intercept=230, color="Urban"), linewidth=1.2) +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  theme_classic() +
  annotate("Text", size=5, label="H1: No change in slope", color="grey20", x=45, y=470)+
  theme(text=element_text(size=13), axis.text.y = element_blank(), axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.position="none")+
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))
#ggsave(h1.plot, file="Desktop/h1.plot.png", height=4.5, width=6)
# H2
h2.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-2.2, intercept=400, linewidth=1.2, color="#009E73") +
  geom_abline(slope=-1.5, intercept=275, linewidth=1.2, color="#CC79A7") +
  geom_abline(slope=-0.75, intercept=150, linewidth=1.2, color="#000000") +
  labs(x="Absolute Latitude", y="Species Richness") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  theme_classic()+
  annotate("Text", size=5, label="H2: Change in slope", color="grey20", x=45, y=470)+
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
 axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position="none")
h2.plot
#ggsave(h2.plot, file="Desktop/h2.plot.png", height=4.5, width=6)

#hypothesis.fig <- ggarrange(h1.plot, h2.plot, labels=c("A","B"), common.legend=TRUE, align="v", legend="right")
#hypothesis.fig <- annotate_figure(hypothesis.fig, bottom=text_grob("Absolute latitude", size = 15, hjust=.8))
#hypothesis.fig
#ggsave(hypothesis.fig, file="committee_meeting_figs/hypothesis_fig.png", height=4, width=8)

##### Seasonal prediction plots

# Winter
winter.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-4.5, intercept=480, color="Natural"), linewidth=1.2) +
  geom_abline(aes(slope=-3.6, intercept=380, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-2.8, intercept=280, color="Urban"), linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
 # annotate("Text", size=5, label="Winter", x=45, y=470, color="grey30") +
  theme_classic() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.title=element_blank(), legend.position="none") +
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  annotate("Text", size=5, label="H3: Winter", color="grey20", x=45, y=470)
winter.plot

#ggsave(winter.plot, file="committee_meeting_figs/pres_winter_hypothesis_fig.png", height=4, width=6)


# Summer

summer.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-2.2, intercept=375, color="Natural"), linewidth=1.2) +
  geom_abline(aes(slope=-1.5, intercept=275, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-0.75, intercept=175, color="Urban"), linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  theme_classic() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.title=element_blank(), legend.position=c(0.85, 0.75)) +
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))+
  annotate("Text", size=5, label="H3: Summer", color="grey20", x=45, y=470)
summer.plot
#ggsave(summer.plot, file="committee_meeting_figs/pres_summer_hypothesis_fig.png", height=4, width=6)


#seasonal.hypothesis.fig <- ggarrange(winter.plot, summer.plot, common.legend=TRUE, align="hv", legend="right")

#seasonal.hypothesis.fig<-annotate_figure(seasonal.hypothesis.fig, bottom=text_grob("Absolute latitude", size = 15, hjust=.9))
#seasonal.hypothesis.fig
#ggsave(seasonal.hypothesis.fig, file="seasonal_hypothesis_fig.png", height=4, width=8)


### Full figure
library(patchwork)

hypothesis.plot <- ggarrange(h1.plot, h2.plot, winter.plot, summer.plot)
require(grid) 
hypothesis.plot2 <- annotate_figure(hypothesis.plot, left = textGrob("Species Richness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                  bottom = textGrob("Absolute Latitude", gp = gpar(cex = 1.3)))
ggsave(hypothesis.plot2, file="hypothesis.plot.png", height=6, width=8)







#### Intermediates
winter.plot2 <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-4.5, intercept=480, color="Natural"), linewidth=1.2) +
#  geom_abline(aes(slope=-3.8, intercept=420, color="Suburban"), linewidth=1.2) +
 # geom_abline(aes(slope=-3, intercept=360, color="Urban"), linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  # annotate("Text", size=5, label="Winter", x=45, y=470, color="grey30") +
  theme_bw() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.position="none") +
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))
ggsave(winter.plot2, file="committee_meeting_figs/pres_winter_hypothesis_fig2.png", height=4, width=6)

winter.plot3 <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-4.5, intercept=480, color="Natural"), linewidth=1.2) +
  #  geom_abline(aes(slope=-3.8, intercept=420, color="Suburban"), linewidth=1.2) +
   geom_abline(aes(slope=-3, intercept=360, color="Urban"), linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  # annotate("Text", size=5, label="Winter", x=45, y=470, color="grey30") +
  theme_bw() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.x=element_blank(), legend.position="none") +
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))
ggsave(winter.plot3, file="committee_meeting_figs/pres_winter_hypothesis_fig3.png", height=4, width=6)



summer.plot2 <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-2, intercept=390, linewidth=1.2, color="#009E73") +
#  geom_abline(slope=-1.6, intercept=330, linewidth=1.2, color="#CC79A7") +
#  geom_abline(slope=-1.2, intercept=270, linewidth=1.2, color="#000000") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  #  annotate("Text", size=5, label="Summer", color="grey30", x=45, y=470) +
  theme_bw() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position="none")
ggsave(summer.plot2, file="committee_meeting_figs/pres_summer_hypothesis_fig2.png", height=4, width=6)


summer.plot3 <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-2, intercept=390, linewidth=1.2, color="#009E73") +
#  geom_abline(slope=-1.6, intercept=330, linewidth=1.2, color="#CC79A7") +
  geom_abline(slope=-1.2, intercept=270, linewidth=1.2, color="#000000") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  #  annotate("Text", size=5, label="Summer", color="grey30", x=45, y=470) +
  theme_bw() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position="none")
ggsave(summer.plot3, file="committee_meeting_figs/pres_summer_hypothesis_fig3.png", height=4, width=6)



############ FULL DATASET ############

#### Get data into place I want it for model
dat <- read.csv("modeling_data.csv")
# turn biome and urban into a factor
dat$BIOME <- as.factor(dat$BIOME)
dat$urban <- as.factor(dat$urban)

dat %>% group_by(urban) %>% summarise(n=n())
summary(dat)
hist(dat$total_SR, breaks=50)

hist(log(dat$total_SR))
hist(sqrt(dat$total_SR), breaks=50) # this looks pretty good
hist(dat$number_checklists) # this is super log normal, used the log in the response variable
dat$abslat <- abs(dat$lat)

dat %>% group_by(BIOME) %>% summarise(n=n())

# make another columbn with only 3 categories
# try model with only 3 categories
dat <- dat %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat$urban2 <- as.factor(dat$urban2)

dat$urban2 <- factor(dat$urban2, levels = c("1", "2", "3"),
                     labels = c("Natural n = 40,490", "Suburban n = 17,623", "Urban n = 12,636"))

# Divide globe into 4 quadrants at the prime meridian and antimeridian
dat$quadrant <- NA

for (i in 1:nrow(dat)){
  if (dat$long[i] < 0 & dat$hemisphere[i] == "northern") { # quadrant 1 is North America
    dat$quadrant[i] <- 1
  }
  else if (dat$long[i] > 0 & dat$hemisphere[i] == "northern") { # quadrant 2 is europe and asia and N Africa
    dat$quadrant[i] <- 2
  }
  else if (dat$long[i] < 0 & dat$hemisphere[i] == "southern") { # quadrant 3 is south america 
    dat$quadrant[i] <- 3 
  }
  else {dat$quadrant[i] <- 4} # quadrant 4 is oceania and southern africa
}
dat$quadrant <- as.factor(dat$quadrant)

dat <- dat %>% filter(!CONTINENT == "Antarctica")

#############
# Plot 1: Mapped coverage of data faceted by urbanization score
world <- ne_countries(scale = "medium", returnclass = "sf")

coverage_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat, aes(x=long, y=lat, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(option="H", na.value = NA, direction=1)+
  labs(x="Longitude", y="Latitude", color="Species \nRichness")+
  theme_bw() +
  facet_wrap(~urban2, ncol=2) +
  theme(legend.position=c(.6,.25))

coverage_plot
ggsave(coverage_plot, file="committee_meeting_figs/coverage_fig.png", height=5, width=8)

###############
# Plot 2: Model results from simple linear model
mod1.full <- lm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                      BIOME + log(number_checklists) + elevation, dat)
#mod1.full <- lm(sqrt(total_SR) ~ abslat * urban2, dat)

predicted <- ggpredict(mod1.full , terms = c("abslat", "urban2")) 

# looks the same whether sqrt included in model or not



results.plot.lm <-
  plot(predicted, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))


results.plot.lm
ggsave(results.plot.lm, file="committee_meeting_figs/results.plot.lm.png", height=5, width=8)
summary(mod1.full)
anova(mod1.full)

library(emmeans)
emmeans(mod1.full, specs="urban2")
emtrends(mod1.full, pairwise ~ urban2, var="abslat")


######################
# Plot 4: Plot results of spatial model with the subset of the data
set.seed(30)
dat.samp <- dat[sample(nrow(dat), 10000), ]

GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtered.tif")
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=10)
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)

dat.samp.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                         BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
saveRDS(dat.samp.sem, "spatialmod.samp.rds")
#dat.samp.sem <- readRDS("spatialmod.samp.rds")
predicted <- predict(dat.samp.sem, interval='confidence')
dat.samp.results <- cbind(dat.samp, predicted)

spat.model.samp <- ggplot(dat.samp.results, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(aes(color=urban2), alpha=0.2)+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_bw()+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'), legend.position="none")+
  labs(x="Absolute Latitude", y="Species Richness")

ggsave(spat.model.samp, file="committee_meeting_figs/spat.model.plot.png", height=5, width=8)






################# SEASONAL DATA ####################
dat.season <- read.csv("season_model_data.csv")
dat.season %>% group_by(season) %>% summarise(n=n())
dat.season$urban<-as.factor(dat.season$urban)
dat.season$BIOME <- as.factor(dat.season$BIOME)

# make another column with just 3 categories
dat.season <- dat.season %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat.season %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat.season$urban2 <- as.factor(dat.season$urban2)
dat.season$abslat <- abs(dat.season$lat) # absolute latitude

# Divide by quartiles
for (i in 1:nrow(dat.season)){
  if (dat.season$long[i] < 0 & dat.season$hemisphere[i] == "northern") { # quadrant 1 is North America
    dat.season$quadrant[i] <- 1
  }
  else if (dat.season$long[i] > 0 & dat.season$hemisphere[i] == "northern") { # quadrant 2 is europe and asia and N Africa
    dat.season$quadrant[i] <- 2
  }
  else if (dat.season$long[i] < 0 & dat.season$hemisphere[i] == "southern") { # quadrant 3 is south america 
    dat.season$quadrant[i] <- 3 
  }
  else {dat.season$quadrant[i] <- 4} # quadrant 4 is oceania and southern africa
}
dat.season$quadrant <- as.factor(dat.season$quadrant)

dat.season <- dat.season %>% filter(!CONTINENT == "Antarctica") # filter out antarctica

dat.season$urban2 <- factor(dat.season$urban2, levels = c("1", "2", "3"),
                     labels = c("Natural", "Suburban", "Urban"))

dat.season$season <- factor(dat.season$season, levels = c("winter", "summer"),
                     labels = c("Winter", "Summer"))

##### Plot 1: Coverage of points for each season
world <- ne_countries(scale = "medium", returnclass = "sf")

season.coverage.plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat.season, aes(x=long, y=lat, color=season), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_manual(values=c("cornflowerblue", "orange"))+
  labs(x="Longitude", y="Latitude") +
  theme_bw()+
  theme(legend.position="none", axis.title = element_blank()) +
  facet_wrap(~season)
ggsave(season.coverage.plot, file="committee_meeting_figs/season.coverage.plot.png", height=5, width=8)


####### Plot 2: Model results from simple linear model
mod1.seas <- lm(sqrt(total_SR) ~ abslat * urban2 * season + quadrant + BIOME + log(number_checklists), dat.season) # including quadruple interaction here

summary(mod1.seas)
anova(mod1.seas)
predicted.season <- ggpredict(mod1.seas, terms = c("abslat", "urban2", "season"))

results.plot.lm.season <-
  plot(predicted.season, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE)) # lines do change if I take out covariates so that's good



ggsave(results.plot.lm.season, file="committee_meeting_figs/results.plot.lm.season.png", height=6, width=12)


predicted <- predict(mod1.seas, interval='confidence')
dat.seas.results <- cbind(dat.season, predicted)

season.model.plot2 <- ggplot(dat.seas.results, aes(x=abs(lat), y=fit, color=urban2))+
  geom_point(aes(color=urban2), alpha=0.2)+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban") +
  facet_wrap(~season) # not sure if this is plotting the right thing

library(emmeans)
# look at model results
emmeans(mod1.seas, specs="urban2", "season")


emtrends(mod1.seas, pairwise ~ urban2, var="abslat", at=c(season="Winter")) # compare slopes in winter
emtrends(mod1.seas, pairwise ~ urban2, var="abslat", at=c(season="Summer")) # compare slopes in summer
emtrends(mod1.seas, pairwise ~ season, var="abslat", at=c(urban2="Natural")) # compare natural slopes in natural

############ Plot 3: Model results from spatial model with subset of data
set.seed(15)
dat.seas.samp <- dat.season[sample(nrow(dat.season), 10000), ]

dat.seas.sf <- st_as_sf(dat.seas.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.seas.nb <- dnearneigh(dat.seas.sf, d1=0, d2=10)
dat.seas.lw <- nb2listw(dat.seas.nb, style = "W", zero.policy = TRUE)


dat.seas.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * season + quadrant + 
                                         BIOME + log(number_checklists), data = dat.seas.samp, listw = dat.seas.lw, zero.policy = TRUE) # run spatial error model

saveRDS(dat.seas.sem, "seasonal.spatialmod.samp.rds")

#dat.seas.sem <- readRDS("seasonal.spatialmod.samp.rds")

predicted <- predict(dat.seas.sem, interval='confidence')
dat.seas.results <- cbind(dat.seas.samp, predicted)

spat.model.seas <- ggplot(dat.seas.results, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(aes(color=urban2), alpha=0.2)+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_bw()+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban") +
  facet_wrap(~season)

ggsave(spat.model.seas, file="committee_meeting_figs/spat.model.plot.season.png", height=6, width=12)


###################### Figs for presentation ################

######## RESULTS OF LINEAR MODEL #############
results.plot.lm2 <-
  plot(predicted, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

ggsave(results.plot.lm2, file="committee_meeting_figs/pres_results.plot.lm.png", height=5, width=9)

# no backtransformation
predicted2 <- ggpredict(mod1.full, terms = c("abslat", "urban2"), back.transform=F) 
?ggeffects::plot


results.plot.lm2 <-
  plot(predicted2, add.data=F, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, ci.style="ribbon") +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

ggsave(results.plot.lm2, file="committee_meeting_figs/pres_results.plot.notrans.png", height=5, width=10)


############## RESULTS OF SEASONAL MODEL #############

results.plot.lm.season <-
  plot(predicted.season, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))


mod1.seas <- lm(sqrt(total_SR) ~ abslat * urban2 * season + quadrant + BIOME + log(number_checklists), dat.season) # including quadruple interaction here

summary(mod1.seas)
anova(mod1.seas)
predicted.season <- ggpredict(mod1.seas, terms = c("abslat", "urban2", "season"))

results.plot.lm.season.pres <-
  plot(predicted.season, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm')) # lines do change if I take out covariates so that's good

ggsave(results.plot.lm.season.pres, file="committee_meeting_figs/pres_results.plot.lm.season.png", height=6, width=12)

# plots with just natural
predicted.season2 <- predicted.season %>% filter(group=="Natural")
results.plot.lm.season.pres2 <-
  plot(predicted.season2, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE)) # lines do change if I take out covariates so that's good

ggsave(results.plot.lm.season.pres2, file="committee_meeting_figs/pres_.lm.season.natural.png", height=6, width=12)

# no back transformation
predicted.season <- ggpredict(mod1.seas, terms = c("abslat", "urban2", "season"), back.transform=F)

results.plot.lm.season <-
  plot(predicted.season, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE)) # lines do change if I take out covariates so that's good

ggsave(results.plot.lm.season, file="committee_meeting_figs/lm.season.notrans.png", height=6, width=12)

### LDG by quadrant 
predicted.quad <- ggpredict(mod1.full , terms = c("abslat", "urban2", "quadrant")) 

results.quad.lm <-
  plot(predicted.quad, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))
?ggpredict

ggsave(results.quad.lm, file="committee_meeting_figs/pres_results.quad.png", height=6, width=10)

#### LDG by continent
mod1.cont <- lm(sqrt(total_SR) ~ abslat * urban2 * CONTINENT + 
                  BIOME + log(number_checklists) + elevation, dat)
#mod1.full <- lm(sqrt(total_SR) ~ abslat * urban2, dat)

predicted.cont <- ggpredict(mod1.cont , terms = c("abslat", "urban2", "CONTINENT")) 

results.cont.lm <-
  plot(predicted.cont, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, ci.style="ribbon",
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'))

ggsave(results.cont.lm, file="committee_meeting_figs/pres_results.cont.png", height=6, width=10)

##### LDG by hemisphere
mod1.hemisphere <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                  BIOME + log(number_checklists) + elevation, dat)

predicted.hemisphere <- ggpredict(mod1.hemisphere, terms = c("abslat", "urban2", "hemisphere")) 

results.hemisphere.lm <-
  plot(predicted.hemisphere, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, ci.style="ribbon",
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'))

ggsave(results.hemisphere.lm, file="committee_meeting_figs/pres_results.hemisphere.png", height=6, width=10)

################ Variogram of spatial autocorrelation
library(nlme)
gls1 <- gls(sqrt(total_SR) ~ abs(lat) * urban * quadrant + 
              BIOME + log(number_checklists) + elevation, dat)


dat.sf <- st_as_sf(dat, coords=c("long", "lat")) 

png(filename="committee_meeting_figs/variogram.png", width=650, height=480, pointsize=40)
plot(gstat::variogram(residuals(gls1, "normalized") ~
                        1, data = dat.sf, cutoff = 100))
dev.off()

