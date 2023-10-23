################ This script is to look into different ways to deal with spatial autocorrelation because it was getting too messy on my main script ####
library(terra)
library(sf)
library(tidyverse)
library(ggeffects)
library(foreach)
library(spdep)
library(spatialreg)
library(elevatr)
library(emmeans)
library(ncf)
library(nlme)

### Load model data
dat <- read.csv("modeling_data.csv")

# Run model
mod1.lm <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   BIOME + log(number_checklists) + elevation, dat)


########### 1. Look at spatial autocorrelation with variogram ########

# rerun using gls 
mod1.gls <- gls(sqrt(total_SR) ~ abs(lat) * urban * hemisphere + 
              BIOME + log(number_checklists) + elevation, dat)


# check for spatial autocorrelation
dat$gls1Resids <- residuals(mod1.gls)
hist(dat$gls1Resids) # looking very normally distributed, that's good

# Tried running a correlogram but it takes up too much memory in R and will not run
#residsI <- spline.correlog(x=dat$long, y=dat$lat,
#      z=dat$gls1Resids, resamp=100, quiet=TRUE) # this takes up too much memory in R because there are so many data points

#plot(residsI,xlim=c(0,50))


#### Variogram

dat.sf <- st_as_sf(dat, coords=c("long", "lat")) 

plot(gstat::variogram(residuals(gls1, "normalized") ~
                        1, data = dat.sf, cutoff = 100))
# big jump at around 60km

# try with directional
plot(gstat::variogram(residuals(gls1, "normalized") ~
                        1, data = dat.sf, cutoff = 100, alpha = c(0, 45, 90, 135)))
# it likely looks like this at 0 because it becomes more similar as you move toward the equator again 

### Moran's I using inverse distance matrix - this is way too computationally intensive and crashes my computer
# make distance matrix
#w <-as.matrix(1/dist(cbind(dat$long, dat$lat))) # this is too much for my computer or the lab computer to handle
#wlist<-mat2listw(w)
#moran.test(dat$gls1Resids, wlist)


dat.nb <- dnearneigh(dat.sf, d1=0, d2=200) # calculate distances
dat.lw <- nb2listw(dat.nb, style = "W", zero.policy = TRUE)
## Run Moran's I
moran.full <- lm.morantest(mod1.lm, dat.lw, zero.policy = T)
beep()
# this is too large for my computer


## Subsample randomly and run Moran's I again
dat.samp <- dat[sample(nrow(dat), 5000), ]
lm.samp <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                BIOME + log(number_checklists) + elevation, dat.samp)
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 
dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=200) # calculate distances
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)
## Run Moran's I
moran.samp <- lm.morantest(lm.samp, dat.samp.lw, zero.policy = T)
moran.samp
# ok, the moran's I value is positive at least, but it is also very 


####### 2. Try to get rid of spatial autocorrelation in the model by fitting a correlation structure in the gls based on the variogram #########

# Fit all possible autocorrelation structures and see which is the best fit
csSpher <- corSpher(form=~lat+long,nugget=TRUE) # sperical
csExp <- corExp(form=~lat+long,nugget=TRUE) # exponential
csGaus <- corGaus(form=~lat+long,nugget=TRUE) # gaussian
csLin <- corLin(form=~lat+long,nugget=TRUE) # linear
csRatio <- corRatio(form=~lat+long,nugget=TRUE) # ratio

# update models
#glsGaus <- update(gls1, correlation=csLin) # gaussian
# This takes too long to run, and crashes my computer


# Try subsampling my data and running the different correlation structures on that
dat.samp <- dat[sample(nrow(dat), 10000), ]
# try to run a correlog

gls1.samp <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + as.factor(BIOME) + log(number_checklists), dat.samp)
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 

plot(gstat::variogram(residuals(gls1.samp, "normalized") ~
                       1, data = dat.samp.sf, cutoff = 100)) # variogram of the sample, looks pretty much the same as the full data
plot(gstat::variogram(residuals(gls1.samp, "normalized") ~
                        1, data = dat.samp.sf, cutoff = 100, alpha = c(0, 45, 90, 135))) # also looks the same as the full one



# Make a correlogram
dat.samp$gls1.samp <- residuals(gls1.samp)
hist(dat.samp$gls1.samp)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$gls1.samp, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,100)) # there is some autocorrelation at small distances

#### Looking at the correlogram for the raw species richness data
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$total_SR, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,100)) # even the raw species richness data does not have a very strong correlation


# calculate Moran's I
w <-as.matrix(1/dist(cbind(dat.samp$long, dat.samp$lat))) # make inverse distance matrix - weights things that are close together higher
wlist<-mat2listw(w) # assign weights based on inverse distances (converts square spatial weights matrix to a weights list object)
# this quickly becomes very large because it is pairwise distances
moran.test(dat.samp$total_SR, wlist)
# it is significant


# update with gaussian
glsGaus <- update(gls1.samp, correlation=csGaus) # gaussian
dat.samp$glsGaus <- residuals(glsGaus)
hist(dat.samp$glsGaus)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$glsGaus, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,100))
beep()
plot(gstat::variogram(residuals(glsGaus, "normalized") ~
                        1, data = dat.samp.sf, cutoff = 100)) # variogram still looks pretty much the same, with the big jump at 60
# still spatially autocorrelated

# update with spherical
glsSpher <- update(gls1.samp, correlation=csSpher)
dat.samp$glsSpher <- residuals(glsSpher)
hist(dat.samp$glsSpher)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$glsSpher, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))

# update with exponential 
glsExp <- update(gls1.samp, correlation=csExp)
dat.samp$glsExp <- residuals(glsExp)
hist(dat.samp$glsExp)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$glsExp, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,100))
plot(gstat::variogram(residuals(glsExp, "normalized") ~
                        1, data = dat.samp.sf, cutoff = 100)) # looks very similar to the gaussian

# update with linear
glsLin <- update(gls1.samp, correlation=csLin)
dat.samp$glsLin <- residuals(glsLin)
hist(dat.samp$glsLin)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$glsLin, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20)) # now the whole thing is autocorrelated, definitely not that

# update with ratio
glsRatio <- update(gls1.samp, correlation=csRatio)
dat.samp$glsRatio <- residuals(glsRatio)
hist(dat.samp$glsRatio)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$glsRatio, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20)) 

# compare AIC
AIC(gls1.samp, glsExp, glsGaus, glsLin, glsSpher, glsRatio)
# glsExp was the best












############ 3. Tried thinning my data to get rid of the spatial autocorrelation ######################


## 3a. I will try to run spThin to a distance where there is no spatial autocorrelation, based on the correlogram of the sample data
## it looks like the autocorrelation ends at about 5 km so I am going to try thinning by that and see what happens

### SpThin only works on samples of the data because it is too computationally intensive so I will take a sample
dat.samp <- dat[sample(nrow(dat), 5000), ]
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 


## Iteratively look at Moran's I to see where the autocorrelation stops
# this takes a really long time
#moran_I <- c()
#for (d in seq(1, 200, 20)) {
#  #foreach  (d = seq(0, 200, 10),
#  #                       .combine = 'c') %dopar% {
#  dat.samp.nb <- dnearneigh(dat.samp.sf, d1 = 0, d2 = d)
#  dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)
#  moran <- moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE)
#  moran_I <- c(moran_I, moran$statistic)
#} # THIS TAKES A REALLY LONG TIME

#beep()
#moran_I <- data.frame(moran = moran_I, 
#                      distance = seq(1, 200, 20))
#
#
#moran <- moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE)
#
#ggplot(moran_I, aes(x = distance, y = moran)) + 
#  geom_point() +
#  geom_line()
#beep()
#

# looking at the variogram when I thin the data. Still looks pretty much the same.
dat.thinned <- thin.algorithm(dat.samp[,c(27:28)], thin.par=100, 
                              rep=1) # this spits out a list of lat long points (list length is # of reps), which I then need to merge back with other data

da.thinned.int <- dat.thinned[[1]] %>% rename(long=Longitude, lat=Latitude)
#
dat.thinned.new <- inner_join(dat.samp, da.thinned.int, by=c("long", "lat"))
dat.thinned.new.sf <- st_as_sf(dat.thinned.new, coords=c("long", "lat")) 
gls.thinned <- gls(sqrt(total_SR) ~ abslat * urban2 * hemisphere + BIOME + log(number_checklists) + elevation, dat.thinned.new)

plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.new.sf, cutoff = 200))
plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.new.sf, cutoff = 100, alpha = c(0, 45, 90, 135), width=1))


## run model with thinned data and check for autocorrelation
lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                  BIOME + log(number_checklists) + elevation, dat.thinned.new)

predicted2 <- ggpredict(lm.thinned, terms = c("abslat", "urban2")) 
## looks the same whether sqrt included in model or not
#
results.plot2 <-
  plot(predicted2, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text = element_text(size=20), legend.spacing.y = unit(1, 'cm'))
results.plot2
#
#
dat.thinned.sf <- st_as_sf(dat.thinned.new, coords=c("long", "lat")) 
dat.thinned.nb <- dnearneigh(dat.thinned.sf, d1=0, d2=1000) # calculate distances
dat.thinned.lw <- nb2listw(dat.thinned.nb, style = "W", zero.policy = TRUE) # turn into weighted list
## supplements a neighbors list with spatial weights for the chosen coding scheme


## Moran's I test
lm.morantest(lm.thinned, dat.thinned.lw, zero.policy = T)
## There is a negative Moran's I value, which indicates that it tends toward dispersion









##### 3b. Trying much faster thinning method
# Create raster grid and overlay and then randomly sample points from the grid
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL)) # sample grid



# assign cell number to each point in my data
vect <- st_as_sf(dat, crs=st_crs(GHSL), coords=c("x","y"))
xy=st_coordinates(vect)
# get cell number that each point is in
dat$cell.subsample<-cellFromXY(sample.grid, xy)

# randomly sample one point within each cell
dat.thinned <- dat %>% group_by(cell.subsample) %>% sample_n(1) 


# run model
lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
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
gls.thinned <- gls(sqrt(total_SR) ~ abslat * urban2 * hemisphere + BIOME + log(number_checklists) + elevation, dat.thinned)
plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.sf, cutoff = 200))

plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.sf, cutoff = 200, alpha = c(0, 45, 90, 135)))

# looking at results for randomly sampled data (not stratified)
dat.samp <- dat[sample(nrow(dat), 3104), ]
lm.samp <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                   BIOME + log(number_checklists) + elevation, dat.samp)

dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 
dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=200) # calculate distances
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE) # turn into weighted list
moran <- lm.morantest(lm.samp, dat.samp.lw, zero.policy = T)
moran # yes this is autocorrelated, so the other one is accurate





## Make map of thinned data
plot.thinned <-  ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat.thinned, aes(x=long, y=lat), color="cornflowerblue", size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude")+
  theme_bw()
plot.thinned







predicted2 <- ggpredict(lm.thinned, terms = c("abslat", "urban2")) 
#
results.plot2 <-
  plot(predicted2, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text = element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

results.plot2
# still getting the same results even when it is super thinned out


summary(lm.thinned)

plot(lm.thinned)



###### Start looping model and store results
thinned.results <- list()
predicted <- list()

for (i in 1:1000){
  dat.thinned <- dat %>% group_by(cell.subsample) %>% sample_n(1) 
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
#  geom_point(aes(x=x, y=conf.high), color="grey40")+
#  geom_point(aes(x=x, y=conf.low), color="grey40")
# natural is definitely different, suburban and urban look the same









################## 4. Trying spatial models in the spatialreg package ###################


# sample data
set.seed(10)
dat.samp <- dat[sample(nrow(dat), 5000), ]
# turn into sf object
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtered.tif")
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   BIOME + log(number_checklists) + elevation, dat.samp) # here is the model that I am working with
summary(mod1.trans)
dat.samp$residuals <- residuals(mod1.trans)
dat.samp$fitted <- fitted(mod1.trans)
#saveRDS(mod1.trans, "50K_samp_lmmod.rds")

dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=200) # make list of nearest neighbors
test <- as.data.frame(card(dat.samp.nb)) # this gives how many neighbors there are. 

class(dat.samp.sf)
# This function identifies neighbours of region points by euclidean distance
# it returns a list of integer vectors giving the region id numbers for neighbors satisfying the distance criteria
class(dat.samp.nb) # this is an nb object

dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)
beep()

# Moran's I test
moran.results <- lm.morantest(mod1.trans, dat.samp.lw, zero.policy = T) # very spatially autocorrelated
moran.results

LMtests.results <- lm.LMtests(mod1.trans, dat.samp.lw, test="all", zero.policy = T) # test for spatial error - very significant
LMtests.results



Inc.lag <- lag.listw(dat.samp.lw, dat.samp$residuals, zero.policy = T)
plot(Inc.lag)
moran.plot(dat.samp$residuals, dat.samp.lw)
# slope of the regression line between spatially lagged values and observed values

# to assess if the slope is significantly diff from 0, can permute values across samples
moran <- moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE) # moran monte carlo
moran
# significant

# look at how it changes with distance of nearest neighbors used
moran_I <- c()


# set up parallelization
n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
# loop d through a sequence ranging from 50 to 2000 ( in parallel)

for (d in seq(1, 200, 10)) {
  #foreach  (d = seq(0, 200, 10),
  #                       .combine = 'c') %dopar% {
  dat.samp.nb <- dnearneigh(dat.samp.sf, d1 = 0, d2 = d)
  dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)
  moran <- moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE)
  moran_I <- c(moran_I, moran$statistic)
} # THIS TAKES A REALLY LONG TIME
beep()
moran_I <- data.frame(moran = moran_I, 
                      distance = seq(0, 200, 10))

ggplot(moran_I, aes(x = distance, y = moran)) + 
  geom_point() +
  geom_line()
beep()



##### Try running models  with sample data
# Spatially lagged X model
dat.samp.slx <- lmSLX(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                        BIOME + log(number_checklists) + elevation, data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
beep()
summary(dat.samp.slx) # R2 is 0.38 (higher than without the lag)
summary(mod1.trans)
summary(impacts(dat.samp.slx, listw=dat.samp.lw), zstats=TRUE) 
#saveRDS(dat.samp.slx, "50K_samp_slxmod.rds")
# interaction effect no longer significant hmmm

# test model
dat.samp$residuals.slx <- residuals(dat.samp.slx)
moran.mc(dat.samp$residuals.slx, dat.samp.lw, nsim = 999, zero.policy = TRUE)
# no more spatial autocorrelation!


# Spatial lag model
dat.samp.slm <- spatialreg::lagsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                       BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
beep()
summary(dat.samp.slm)

# Spatial error model
dat.samp.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                         BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
summary(dat.samp.sem)
beep()


AIC(dat.samp.slx, dat.samp.sem, dat.samp.slm)
# the slx model is best

# run moran test
dat.samp$residuals.slx <- residuals(dat.samp.slx)
dat.samp$residuals.slm <- residuals(dat.samp.slm)
dat.samp$residuals.sem <- residuals(dat.samp.sem)

moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE) # significant p = 0.001
moran.mc(dat.samp$residuals.slx, dat.samp.lw, nsim = 999, zero.policy = TRUE) # still autocorrelated
moran.mc(dat.samp$residuals.slm, dat.samp.lw, nsim = 999, zero.policy = TRUE) # still autocorrelated
moran.mc(dat.samp$residuals.sem, dat.samp.lw, nsim = 999, zero.policy = TRUE) # not autocorrelated! Spatial error model could be the answer!
beep()
# the spatial error model is the best!


####### Try different distances
set.seed(15)
dat.samp2 <- dat[sample(nrow(dat), 5000), ]

# run regular model
dat.samp2.lm <- lm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                     BIOME + log(number_checklists), data = dat.samp2)

# turn into sf object
dat.samp2.sf <- st_as_sf(dat.samp2, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.samp2.nb <- dnearneigh(dat.samp2.sf, d1=0, d2=5)
dat.samp2.lw <- nb2listw(dat.samp2.nb, style = "W", zero.policy = TRUE)

dat.samp2.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                          BIOME + log(number_checklists), data = dat.samp2, listw = dat.samp2.lw, zero.policy = TRUE)

dat.samp2$residuals.lm <- residuals(dat.samp2.lm)
moran.mc(dat.samp2$residuals.lm, dat.samp2.lw, nsim = 999, zero.policy = TRUE)


dat.samp2$residuals.sem <- residuals(dat.samp2.sem)
moran.mc(dat.samp2$residuals.sem, dat.samp2.lw, nsim = 999, zero.policy = TRUE) # 100 km still gets rid of autocorrelation
# still works with a distance of 50 and 10 and 1!
beep()

plot(ggeffects::ggpredict(dat.samp2.sem, terms=c("abslat"), listw=dat.samp2.lw, facets = TRUE))

summary(dat.samp2.sem)
# don't use impacts function on error models

plot(predict(dat.samp2.sem, interval='confidence')[1:25])
predicted <- predict(dat.samp2.sem, interval='confidence')
dat.samp2.test <- cbind(dat.samp2, predicted)

ggplot(dat.samp2.test, aes(x=abs(lat), y=fit^2, color=urban2))+
  geom_point(alpha=0.1)+
  geom_smooth(method="lm") # plot model fits yay!
# the model results still look the same

################ Try running model on full data
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtered.tif")
dat.sf <- st_as_sf(dat, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.nb <- dnearneigh(dat.sf, d1=0, d2=1)
dat.lw <- nb2listw(dat.nb, style = "W", zero.policy = TRUE)

dat.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                    BIOME + log(number_checklists), data = dat, listw = dat.lw, zero.policy = TRUE)
# can't run, sample size too large
dat.samp2$residuals.lm <- residuals(dat.samp2.lm)
moran.mc(dat.samp2$residuals.lm, dat.samp2.lw, nsim = 999, zero.policy = TRUE)






