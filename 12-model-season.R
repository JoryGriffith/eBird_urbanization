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

# Plot relationship 

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



################################
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


################## Loop and store models

thinned.results.season <- list()
predicted.season <- list()
means.season <- list()
emmeans.slopes.sum <- list()
emmeans.slopes.wint <- list()
ggeffects.slopes.season <- list()
set.seed(20)
for (i in 1:1000){
  dat.thinned <- dat %>% group_by(cell.subsample, season) %>% sample_n(1) 
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere +
                     BIOME + log(number_checklists) + elevation, dat.thinned)
  
  thinned.results.season[[i]] <- summary(lm.thinned) # store summary data
  predicted.season[[i]] <- ggemmeans(lm.thinned, terms = c("abslat", "urban2", "season"))# store predictions
  means.season[[i]] <- emmeans(lm.thinned, specs="urban2", by="season")
  emmeans.slopes.sum[[i]] <- emtrends(lm.thinned, pairwise ~ urban2, var="abslat", at=c(season="Summer")) # so I can see differences in slopes for each model (using emmeans)
  emmeans.slopes.wint[[i]] <- emtrends(lm.thinned, pairwise ~ urban2, var="abslat", at=c(season="Winter"))
  ggeffects.slopes.season[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2", "season"), test=NULL) # see differences in slopes (using ggeffects)
}

# loop and store models

## Save predicted values as a csv
predicted.season.df <- bind_rows(predicted.season)
write.csv(predicted.season.df, "thinned.seasonal.results.csv")

# 1) Look at mean species richness in each urbanization level and season 
seasonal_means_df <- list()
for (i in 1:1000){
  seasonal_means_df[[i]] <- as.data.frame(season.means[[i]])
}

seasonal_means_df <- bind_rows(seasonal_means_df)
seasonal_means_summary <- seasonal_means_df %>% group_by(season, urban2) %>% summarise(mean=mean(emmean), max.upper=mean(upper.CL), min.lower=min(lower.CL))
seasonal_means_summary

# 2) Look at which slopes are different from one another
contrast_sum_df <- list()
for (i in 1:1000){
  contrast_sum_df[[i]] <- as.data.frame(emmeans.slopes.sum[[i]]$contrasts)
}
contrast_sum_df <- bind_rows(contrast_sum_df)
# calculate proportion where each contrast is significant
contrast_sum_df %>% filter(p.value<0.05) %>% group_by(contrast) %>% count()
# suburban and urban are not significantly different but natural is different from urban and suburban (P-value is less than 0.05 in all models)

## Look at the slope of each line
contrast_wint_df <- list()
for (i in 1:1000){
  contrast_wint_df[[i]] <- as.data.frame(emmeans.slopes.wint[[i]]$contrasts)
}
contrast_wint_df <- bind_rows(contrast_wint_df)
contrast_wint_df %>% filter(p.value<0.05) %>% group_by(contrast) %>% count()



### 3) Look at slopes of the line
## Do with both emmeans and ggeffects because they give different results

## Emmeans
slopes_sum_df <- list()
for (i in 1:1000){
  slopes_sum_df[[i]] <- as.data.frame(emmeans.slopes.sum[[i]]$emtrends)
}
slopes_sum_df <- bind_rows(slopes_sum_df)

## Look at the slope of each line
slopes_wint_df <- list()
for (i in 1:1000){
  slopes_wint_df[[i]] <- as.data.frame(emmeans.slopes.wint[[i]]$emtrends)
}
slopes_wint_df <- bind_rows(slopes_wint_df)

# look at slopes in urban areas
slopes_sum <- slopes_sum_df %>% group_by(urban2) %>% summarise(mean=mean(abslat.trend), conf.high = max(upper.CL), conf.low=min(lower.CL))
slopes_sum
slopes_wint <- slopes_wint_df %>% group_by(urban2) %>% summarise(mean=mean(abslat.trend), 
                                                                 conf.high = max(upper.CL), conf.low=min(lower.CL))
slopes_wint
# slopes are significantly negative in both winter and summer, but barely in urban areas in the winter

## Ggeffects
ggeffects.slopes.season.df <- bind_rows(ggeffects.slopes.season)
# this uses the marginaleffects package
# the slopes are negative here for urban and suburban

ggslopes_sum <- ggeffects.slopes.season.df %>% group_by(urban2, season) %>% 
  summarise(mean=mean(Slope), conf.high = max(conf.high), conf.low=min(conf.low))
ggslopes_sum












