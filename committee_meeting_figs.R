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
  facet_wrap(~urban2, nrow=3)

coverage_plot
ggsave(coverage_plot, file="committee_meeting_figs/coverage_fig.png")

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
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

results.plot.lm
ggsave(results.plot.lm, file="committee_meeting_figs/results.plot.lm.png", height=5, width=9)


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
ggsave(results.plot.lm, file="committee_meeting_figs/results.plot.lm.quad.png", height=8, width=9)

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

predicted <- predict(dat.samp.sem, interval='confidence')
dat.samp.results <- cbind(dat.samp, predicted)

ggplot(dat.samp.results, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(aes(color=urban2), alpha=0.1)+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")



















