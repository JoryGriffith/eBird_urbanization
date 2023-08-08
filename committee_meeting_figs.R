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
  theme_bw() +
  theme(text=element_text(size=13), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.x=element_blank())+
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
  theme_bw()+
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
 axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())
#ggsave(h2.plot, file="Desktop/h2.plot.png", height=4.5, width=6)

hypothesis.fig <- ggarrange(h1.plot, h2.plot, labels=c("A","B"), common.legend=TRUE, align="v", legend="right")
hypothesis.fig <- annotate_figure(hypothesis.fig, bottom=text_grob("Absolute latitude", size = 15, hjust=.8))
hypothesis.fig
ggsave(hypothesis.fig, file="committee_meeting_figs/hypothesis_fig.png", height=4, width=8)
?ggarrange

##### Seasonal prediction plots

# Winter
winter.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(aes(slope=-4.5, intercept=480, color="Natural"), linewidth=1.2) +
  geom_abline(aes(slope=-4, intercept=420, color="Suburban"), linewidth=1.2) +
  geom_abline(aes(slope=-3.5, intercept=360, color="Urban"), linewidth=1.2) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  annotate("Text", size=5, label="Winter", x=45, y=470, color="grey30") +
  theme_bw() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(name="Urbanization", breaks=c('Natural', 'Suburban', 'Urban'),
                     values=c('Natural'='#009E73', 'Suburban'='#CC79A7', 'Urban'='#000000'))


# Summer

summer.plot <- ggplot(dat, aes(x=x, y=y))+
  geom_point(alpha=0)+
  geom_abline(slope=-2, intercept=390, linewidth=1.2, color="#009E73") +
  geom_abline(slope=-1.2, intercept=310, linewidth=1.2, color="#CC79A7") +
  geom_abline(slope=-0.4, intercept=230, linewidth=1.2, color="#000000") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,500), expand=c(0,0))+
  labs(x="Absolute Latitude", y="Species Richness") +
  annotate("Text", size=5, label="Summer", color="grey30", x=45, y=470) +
  theme_bw() +
  theme(text=element_text(size=15), axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())

seasonal.hypothesis.fig <- ggarrange(winter.plot, summer.plot, common.legend=TRUE, align="hv", legend="right")

seasonal.hypothesis.fig<-annotate_figure(seasonal.hypothesis.fig, bottom=text_grob("Absolute latitude", size = 15, hjust=.9))
seasonal.hypothesis.fig
ggsave(seasonal.hypothesis.fig, file="committee_meeting_figs/seasonal_hypothesis_fig.png", height=4, width=8)

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

library(emmeans)
emmeans(mod1.full, specs="urban2")
emtrends(mod1.full, specs="urban2", var="abslat") # need to back transform this to make it meaningful

11.69^2-10.69^2
10.69^2
11.69^2-9.54^2
###################
# Plot 3: Plot results of model by quadrant
predicted2 <- ggpredict(mod1.full , terms = c("abslat", "urban2", "quadrant")) 
# looks the same whether sqrt included in model or not


results.plot.lm.quad <-
  plot(predicted2, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

results.plot.lm.quad
ggsave(results.plot.lm, file="committee_meeting_figs/results.plot.lm.quad.png", height=8, width=10)

######################
# Plot 4: Plot results of spatial model with the subset of the data
set.seed(30)
dat.samp <- dat[sample(nrow(dat), 10000), ]

GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtered.tif")
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=200)
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)

dat.samp.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                         BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
saveRDS(dat.samp.sem, "spatialmod.samp.rds")
dat.samp.sem <- readRDS("spatialmod.samp.rds")
predicted <- predict(dat.samp.sem, interval='confidence')
dat.samp.results <- cbind(dat.samp, predicted)

spat.model.samp <- ggplot(dat.samp.results, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(aes(color=urban2), alpha=0.2)+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")

ggsave(spat.model.samp, file="committee_meeting_figs/spat.model.plot.png", height=8, width=10)






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
mod1.seas <- lm(sqrt(total_SR) ~ abslat * urban2 * season * quadrant +
             abslat:hemisphere + BIOME + log(number_checklists), dat.season) # including quadruple interaction here

predicted.season <- ggpredict(mod1.seas, terms = c("abslat", "urban2", "season"))

results.plot.lm.season <-
  plot(predicted.season, add.data=TRUE, dot.size=0.5, alpha=.6, dot.alpha=0.4, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000"), show.legend=FALSE, facet=TRUE) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE)) 

ggsave(results.plot.lm.season, file="committee_meeting_figs/results.plot.lm.season.png", height=8, width=12)


############ Plot 3: Model results from spatial model with subset of data
set.seed(15)
dat.seas.samp <- dat.season[sample(nrow(dat.season), 10000), ]

dat.seas.sf <- st_as_sf(dat.seas.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.seas.nb <- dnearneigh(dat.seas.sf, d1=0, d2=5)
dat.seas.lw <- nb2listw(dat.seas.nb, style = "W", zero.policy = TRUE)


dat.seas.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * season * quadrant + 
                                         BIOME + log(number_checklists), data = dat.seas.samp, listw = dat.seas.lw, zero.policy = TRUE) # run spatial error model

saveRDS(dat.seas.sem, "seasonal.spatialmod.samp.rds")

predicted <- predict(dat.seas.sem, interval='confidence')
dat.seas.results <- cbind(dat.seas.samp, predicted)

spat.model.samp <- ggplot(dat.seas.results, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(aes(color=urban2), alpha=0.2)+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban") +
  facet_wrap(~season)

ggsave(spat.model.samp, file="committee_meeting_figs/spat.model.plot.png", height=8, width=14)





