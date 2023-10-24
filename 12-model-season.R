########### This is the script for modelling the seasonal differences in bird diversity ###################
library(tidyverse)
library(terra)
library(sp)
library(nlme)
library(ncf)
library(rnaturalearth)
library(rnaturalearthdata)
library(spdep)
library(sf)
library(spatialreg)
library(mgcv)
library(ggeffects)

world <- ne_countries(scale = "medium", returnclass = "sf")
# Load data
dat <- read.csv("season_modeling_data.csv")

# assign hemisphere
dat$hemisphere <- "northern"
dat$hemisphere[dat$lat<0]<-"southern"

dat$urban<-as.factor(dat$urban)
dat$BIOME <- as.factor(dat$BIOME)

# make another column with just 3 categories
dat <- dat %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat$urban2 <- as.factor(dat$urban2)
dat$abslat <- abs(dat$lat) # absolute latitude

# Divide by quartiles
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

#dat <- dat %>% filter(!CONTINENT == "Antarctica") # filter out antarctica


urb <- dat %>% filter(urban2==3)
range(urb$lat) # the highest latitude is 64.15 and the lowest is -54.84

suburb <- dat %>% filter(urban2==2)
range(suburb$lat) # highest latitude is 65.76 and lowest latitude is -53.11

nat <- dat %>% filter(urban2==1)
range(nat$lat) # highest latitude is 65.76 and lowest latitude is -53.11
# 42266 - 42243
dat <- dat %>% filter(lat <= 66 & lat >= -55) # cut off data at the highest latitude range

nrow(dat %>% filter(season=="winter"))

dat %>% group_by(season, urban2) %>% count()
# Plot relationship 

LDG <- ggplot(dat, aes(x=abs(lat), y=total_SR, color=urban2)) +
  geom_point(alpha=0.1, shape=1) +
  geom_smooth(method="lm") +
  labs(x="Absolute Latitude", y="Species Richness")+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw() +
  facet_wrap(~season) # has the expected relationship!

# Plot data coverage
dat.wint  <- dat %>% filter(season=="winter") 
plot.wint <-  ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat.wint, aes(x=long, y=lat), color="cornflowerblue", size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude")+
  theme_bw()
ggsave(plot.wint, file="winter.coverage.png", height=8, width=6)

dat.sum  <- dat %>% filter(season=="summer") 
plot.sum <-  ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat.sum, aes(x=long, y=lat), color="orange", size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude")+
  theme_bw()
ggsave(plot.sum, file="summer.coverage.png", height=8, width=6)

dat$season <- factor(dat$season, levels = c("winter", "summer"),
                  labels = c("Winter", "Summer"))

dat$urban2 <- factor(dat$urban2, levels = c("1", "2", "3"),
                     labels = c("Natural", "Suburban", "Urban"))
# model
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + BIOME + log(number_checklists) + elevation, dat)

#mod2 <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + BIOME + log(number_checklists) + elevation, dat)
AIC(mod1, mod2) # better with the interaction
summary(mod1)
anova(mod1) # triple interaction is significant

# look at model residuals
plot(mod1)
# looking pretty good

# run a gls
gls1 <- gls(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere +
              abslat:hemisphere + BIOME + log(number_checklists), dat)


# look at autocorrelation
datSPDF <- dat
coordinates(datSPDF) <- c("long","lat")
plot(gstat::variogram(residuals(gls1) ~
                        1, data = datSPDF, cutoff = 100))
# looks pretty similar

# try directional
plot(gstat::variogram(residuals(gls1, "normalized") ~
                        1, data = datSPDF, cutoff = 100, alpha = c(0, 45, 90, 135)))
# definitely looks autocorrelated

# some plots
library(jtools)
effect_plot(mod1, abslat)
effect_plot(mod1, hemisphere) # higher richness is N
effect_plot(mod1, BIOME)
effect_plot(mod1, number_checklists)
effect_plot(mod1, season) # higher richness in summer, good

library(interactions)
interact_plot(mod1, abslat, urban2, interval=TRUE) # whoa
interact_plot(mod1, abslat, season, interval=TRUE)
interact_plot(mod1, abslat, hemisphere, interval=TRUE) # interesting, it is steeper in the southern hemisphere
# steeper gradient in winter as predicted

mod.test <- lm(sqrt(total_SR) ~ abslat * urban2 * season, dat)
preds <- (predict(mod.test))^2 # predict and back transform
# plot model predictions
ggplot(dat, aes(x=abslat, y=sqrt(total_SR), color=urban2)) + 
  facet_wrap(~season)+
  geom_line(data = cbind(dat, preds), aes(y = preds)) # it is looking fucked up because there are so many other variables that i am not including in the plot

summary(mod1)
anova(mod1)
#  geom_point(alpha=0.1)
library(emmeans)
library(ggeffects)

####### Results
emmeans(mod1, specs=c("urban2", "season")) 
emmeans(mod1, specs=c("abslat", "urban2", "season")) 
lstrends(mod1, pairwise ~ season, var="abslat", by="urban2")

lstrends(mod1, pairwise ~ urban2, var="abslat", at=c(season="Winter")) # compare slopes in winter
lstrends(mod1, pairwise ~ urban2, var="abslat", at=c(season="Summer")) # compare slopes in summer
# still significantly positive!



predicted <- ggpredict(mod1, terms = c("abslat", "urban2", "season")) # looks the same whether sqrt included in model or not

seasonal.results.plot<-
  plot(predicted, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.text=element_text(size=22), legend.spacing.y = unit(1, 'cm'), legend.position="none")+
  guides(fill = guide_legend(byrow = TRUE))
seasonal.results.plot
ggsave(seasonal.results.plot, file="seasonal.results.plot.png", height=6, width=12)
ggsave(seasonal.results.plot, file="seasonal.results.plot.conclusion.png", height=6, width=7)

predicted2 <- predicted %>% filter(group=="Natural")
seasonal.results.plot2<-
  plot(predicted2, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#009E73")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.text=element_text(size=22), legend.spacing.y = unit(1, 'cm'), legend.position="none")+
  guides(fill = guide_legend(byrow = TRUE))
seasonal.results.plot2
ggsave(seasonal.results.plot2, file="seasonal.results2.plot.png", height=6, width=12)
ggsave(seasonal.results.plot2, file="seasonal.results2.plot.conclusion.png", height=6, width=7) # less wide for conclusion

predicted3 <- predicted %>% filter(group%in%c("Natural", "Urban"))
seasonal.results.plot3<-
  plot(predicted3, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.text=element_text(size=22), legend.spacing.y = unit(1, 'cm'), legend.position="none")
seasonal.results.plot3
ggsave(seasonal.results.plot3, file="seasonal.results.plot3.png", height=6, width=12)
ggsave(seasonal.results.plot3, file="seasonal.results3.plot.conclusion.png", height=6, width=7)
# it back transforms automatically but does not transform the errors, might need to fix later



# subsample and run model to look at correlation
dat.samp <- dat[sample(nrow(dat), 1000), ]

gls1.samp <- gls(sqrt(total_SR) ~ abslat * urban2 * season + hemisphere +
              abslat:hemisphere + BIOME + log(number_checklists), dat.samp)
summary(gls1.samp)
anova(gls1.samp)
# look at autocorrelation
dat.samp$gls1.samp <- residuals(gls1.samp)
hist(dat.samp$gls1.samp)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$gls1.samp, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))

# try to add a correlation structure
gls1.samp.cor <- update(gls1, correlation=corSpatial(form = ~ lat + long | season , nugget=TRUE)) 
beep()
# need to add a grouping factor for season - some overlapping points in different seasons
?gls
dat.samp$gls1.samp.cor <- residuals(gls1.samp.cor)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$gls1.samp.cor, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
AIC(gls1.samp, gls1.samp.cor)



########################
# Try new spatial models

# try to run in parallel
library(parallel)
nc <- detectCores(logical=FALSE)
cl <- makeCluster(nc)
set.ClusterOption(cl)
get.ClusterOption()

dat.samp <- dat[sample(nrow(dat), 37057), ]

GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtered.tif")
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

lm1.samp <- lm(sqrt(total_SR) ~ abslat * urban2 * season + quadrant 
                 + BIOME + log(number_checklists) + elevation, dat.samp)

dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=5)
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)

lm.morantest(lm1.samp, dat.samp.lw, zero.policy = T) # test for autocorrelation - signficant
lm.LMtests(lm1.samp, dat.samp.lw, test="all", zero.policy = T) # test for spatial error - very significant

dat.samp.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * season + quadrant + 
                                         BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE) # run spatial error model

dat.samp$residuals.sem <- residuals(dat.samp.sem)
moran.mc(dat.samp$residuals.sem, dat.samp.lw, nsim = 999, zero.policy = TRUE) # autocorrelation test - not sig!

# plot results
predicted <- predict(dat.samp.sem, interval='confidence')
dat.samp.test <- cbind(dat.samp, predicted)

ggplot(dat.samp.test, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(alpha=0.1)+
  geom_smooth(method="lm") +
  facet_wrap(~season)
# similar results yay!


### Try to run with full data
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtered.tif")
dat.sf <- st_as_sf(dat, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.nb <- dnearneigh(dat.sf, d1=0, d2=5)
dat.lw <- nb2listw(dat.nb, style = "W", zero.policy = TRUE)

dat.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * season + quadrant + 
                                         BIOME + log(number_checklists), data = dat.samp, listw = dat.lw, zero.policy = TRUE) # run spatial error model



########################
### Try non-linear models

## Northern hemisphere winter
dat.nw <- dat %>% filter(hemisphere=="northern" & season=="Winter") 

mod.gam1 <- gam(total_SR ~ s(lat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data=dat.nw)
plot(ggeffects::ggpredict(mod.gam1, terms=c("lat", "urban2"), facets = TRUE))


#### N hemsiphere summer
dat.ns <- dat %>% filter(hemisphere=="northern" & season=="Summer") 

mod.gam2 <- gam(total_SR ~ s(lat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data=dat.ns)
plot(ggeffects::ggpredict(mod.gam2, terms=c("lat", "urban2"), facets = TRUE, add.data=TRUE))


### South hemisphere winter
dat.sw <- dat %>% filter(hemisphere=="southern" & season=="Winter") 

mod.gam3 <- gam(total_SR ~ s(lat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data=dat.sw)
plot(ggeffects::ggpredict(mod.gam3, terms=c("lat", "urban2"), facets = TRUE))

#### S hemisphere summer
dat.sn <- dat %>% filter(hemisphere=="southern" & season=="Summer") 
mod.gam4 <- gam(total_SR ~ s(lat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data=dat.sn)
plot(ggeffects::ggpredict(mod.gam4, terms=c("lat", "urban2"), facets = TRUE), add.data=TRUE)
# pretty much the same pattern



############################
## Iterative thinned models to get rid of spatial autocorrelation
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL)) # sample grid



# assign cell number to each point in my data
vect <- st_as_sf(dat, crs=st_crs(GHSL), coords=c("x","y"))
xy=st_coordinates(vect)
# get cell number that each point is in
dat$cell.subsample<-cellFromXY(sample.grid, xy)


dat.thinned <- dat %>% group_by(cell.subsample, season) %>% sample_n(1) 


# run model
lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere +
                   BIOME + log(number_checklists) + elevation, dat.thinned)
dat.thinned$residuals <- residuals(lm.thinned)

dat.thinned.sf <- st_as_sf(dat.thinned, coords=c("long", "lat")) 
dat.thinned.nb <- dnearneigh(dat.thinned.sf, d1=0, d2=200) # calculate distances
dat.thinned.lw <- nb2listw(dat.thinned.nb, style = "W", zero.policy = TRUE) # turn into weighted list
# supplements a neighbors list with spatial weights for the chosen coding scheme
# trying with B instead of W


# Moran's I test
moran <- lm.morantest(lm.thinned, dat.thinned.lw, zero.policy = T)

moran
# thinned by a spatial grid of 20, the p value is 1, the observed Moran's I is very small. 14269 observations.
# Tried by a spatial grid of 10, still not significantly autocorrelated. This is what I will use

#gls.thinned <- gls(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
#                  BIOME + log(number_checklists) + elevation, dat.thinned)
gls.thinned <- gls(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + BIOME + log(number_checklists) + elevation, dat.thinned)
plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.sf, cutoff = 200))

plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.sf, cutoff = 200, alpha = c(0, 45, 90, 135)))

## might need to thin within season

thinned.results <- list()
predicted <- list()

for (i in 1:1000){
  dat.thinned <- dat %>% group_by(cell.subsample, season) %>% sample_n(1) # try subsampling within season
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                     BIOME + log(number_checklists) + elevation, dat.thinned)
  # dat.thinned.sf <- st_as_sf(dat.thinned, coords=c("long", "lat")) 
  #dat.thinned.nb <- dnearneigh(dat.thinned.sf, d1=0, d2=200) # calculate distances
  #  dat.thinned.lw <- nb2listw(dat.thinned.nb, style = "W", zero.policy = TRUE) # turn into weighted list
  # moran <- lm.morantest(lm.thinned, dat.thinned.lw, zero.policy = T)
  thinned.results[[i]] <- summary(lm.thinned)
  predicted[[i]] <- ggpredict(lm.thinned, terms = c("abslat", "urban2")) 
}

predicted_df <- bind_rows(predicted)

# plot each predicted value as a point and the confidence intervals as lines
predicted_df <- predicted_df %>% group_by(x, group) %>% mutate(max.conf.high = max(conf.high), min.conf.low = min(conf.low))

ggplot(predicted_df, aes(x=x, y=predicted, color=group)) +
  geom_point()+
  geom_smooth(method="lm") +
  geom_errorbar(aes(ymin=min.conf.low, ymax=max.conf.high))




