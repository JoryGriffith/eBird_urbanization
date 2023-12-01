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
library(marginaleffects)

### Load model data
dat <- read.csv("modeling_data.csv")

# Run model
mod1.lm <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat)


########### 1. Look at spatial autocorrelation with variogram ########

# rerun using gls 
mod1.gls <- gls(sqrt(total_SR) ~ abs(lat) * urban * hemisphere + 
              precip + log(number_checklists) + elevation, dat)


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
                precip + log(number_checklists) + elevation, dat.samp)
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
                   abs(lat):CONTINENT + as.factor(precip) + log(number_checklists), dat.samp)
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
gls.thinned <- gls(sqrt(total_SR) ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, dat.thinned.new)

plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.new.sf, cutoff = 200))
plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.new.sf, cutoff = 100, alpha = c(0, 45, 90, 135), width=1))


## run model with thinned data and check for autocorrelation
lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                  precip + log(number_checklists) + elevation, dat.thinned.new)

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
                   precip + log(number_checklists) + elevation, dat.thinned)
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
#                  precip + log(number_checklists) + elevation, dat.thinned)
gls.thinned <- gls(sqrt(total_SR) ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, dat.thinned)
plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.sf, cutoff = 200))

plot(gstat::variogram(residuals(gls.thinned, "normalized") ~
                        1, data = dat.thinned.sf, cutoff = 200, alpha = c(0, 45, 90, 135)))

# looking at results for randomly sampled data (not stratified)
dat.samp <- dat[sample(nrow(dat), 3104), ]
lm.samp <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                   precip + log(number_checklists) + elevation, dat.samp)

dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 
dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=200) # calculate distances
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE) # turn into weighted list
moran <- lm.morantest(lm.samp, dat.samp.lw, zero.policy = T)
moran # yes this is autocorrelated, so the other one is accurate







###### Start looping model and store results
square <- function(x){
  x^2
} # make function to square
thinned.results <- list()
predicted <- list()
means <- list()
emmeans.slopes <- list()
ggeffects.slopes <- list()
ggeffects.slopes.contrast <- list()
ggeffects.slopes.contrast.hemisphere <- list()
set.seed(30)

for (i in 1:1000){
  dat.thinned <- dat %>% group_by(cell.subsample) %>% sample_n(1) 
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned)
  # dat.thinned.sf <- st_as_sf(dat.thinned, coords=c("long", "lat")) 
  #dat.thinned.nb <- dnearneigh(dat.thinned.sf, d1=0, d2=200) # calculate distances
  #  dat.thinned.lw <- nb2listw(dat.thinned.nb, style = "W", zero.policy = TRUE) # turn into weighted list
  # moran <- lm.morantest(lm.thinned, dat.thinned.lw, zero.policy = T)
  thinned.results[[i]] <- anova(lm.thinned) # store summary data
  predicted[[i]] <- avg_predictions(lm.thinned, by=c("abslat", "urban2"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"))) # store predictions for plotting
  means[[i]] <- marginal_means(lm.thinned, variables=c("abslat", "urban2"), transform=square)
  #  means[[i]] <- emmeans(lm.thinned, specs="urban2")
  emmeans.slopes[[i]] <- emtrends(lm.thinned, pairwise ~ urban2, var="abslat") # so I can see differences in slopes for each model (using emmeans)
  ggeffects.slopes[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"), test = NULL) # see differences in slopes (using ggeffects)
  ggeffects.slopes.contrast[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"))
  ggeffects.slopes.contrast.hemisphere[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2", "hemisphere"))
  }

predicted_df <- bind_rows(predicted)
write.csv(predicted_df, "thinned.results.csv")

saveRDS(thinned.results, file="thinned_results/thinned_anovas.rds")
#test <- readRDS("thinned_results/thinned_anovas.rds")

# 1) Look at mean species richness in each urbanization level
means_df <- bind_rows(means)
write.csv(means_df, file="thinned_results/thinned_means.csv")
means_df <- read.csv("thinned_results/thinned_means.csv")
means_summary <- means_df %>% group_by(urban2) %>% summarise(mean=mean(estimate), max.upper=mean(conf.high), min.lower=min(conf.low))
means_summary
129-106
129-92.6


# 2) Look at which slopes are different from one another
slopes_df <- list()
for (i in 1:1000){
  slopes_df[[i]] <- as.data.frame(emmeans.slopes[[i]]$contrasts)
}

slopes_df2 <- bind_rows(slopes_df)
contrast <- slopes_df2 %>% filter(contrast=="Suburban - Urban") %>% filter(p.value<0.05) # 87.9 of the urban-surburban contrasts 
# have p-values that are significantly different from one another
contrast <- slopes_df2 %>% filter(contrast=="Natural - Suburban") %>% filter(p.value<0.05) # 100% of models are significant
contrast <- slopes_df2 %>% filter(contrast=="Natural - Urban") %>% filter(p.value<0.05) # 100% of models are significant



# 3) Look at the slopes of the line
## See if the ggeffects slopes are different than the lsmeans
# lsmeans slopes
emmeans.slopes.df <- list()
for (i in 1:1000){
  emmeans.slopes.df[[i]] <- as.data.frame(emmeans.slopes[[i]]$emtrends)
}
emmeans.slopes.df <- bind_rows(emmeans.slopes.df)
emmeans_slopes <- emmeans.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(abslat.trend), conf.high = max(upper.CL), conf.low=min(lower.CL))
emmeans_slopes

# ggeffects slopes
ggeffects.slopes.df <- bind_rows(ggeffects.slopes)
write.csv(ggeffects.slopes.df, file="thinned_results/thinned_slopes.csv")
ggeffects.slopes.df <- read.csv("thinned_results/thinned_slopes.csv")
ggslopes <- ggeffects.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), conf.low=min(conf.low))
ggslopes


ggeffects.contrast_df <- bind_rows(ggeffects.slopes.contrast)
write.csv(ggeffects.contrast_df, file="thinned_results/thinned_contrasts.csv")
contrast <- ggeffects.contrast_df %>% filter(urban2=="Urban-Suburban") %>% filter(p.value<0.05) # 918, just shy of significant
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Suburban") %>% filter(p.value<0.05) #1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Urban") %>% filter(p.value<0.05) # 1000

# look at contrasts for latitude
ggeffects.contrast.hemisphere_df <- bind_rows(ggeffects.slopes.contrast.hemisphere)
write.csv(ggeffects.contrast.hemisphere_df, file="thinned_results/thinned_constrasts_hemisphere.csv")

contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="Natural-Natural", hemisphere=="northern-southern") %>% filter(p.value<0.05)
# steeper in N hemsiphere
contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="Urban-Urban", hemisphere=="northern-southern") %>% filter(p.value<0.05)
contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="Suburban-Suburban", hemisphere=="northern-southern") %>% filter(p.value<0.05)

contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="Urban-Suburban", hemisphere=="southern-southern") %>% filter(p.value<0.05)
# 718
contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="Urban-Suburban", hemisphere=="northern-northern") %>% filter(p.value<0.05)

## Run full model and plot model fits from that with these confidence intervals
#predicted.full <- ggpredict(mod1.lm, terms = c("abslat", "urban2")) 
#predicted.full$
#
#results.plot <-
#  plot(predicted.full, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
#       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
#  theme_bw()+
#  labs(x="Absolute Latitude", y="Species Richness", color="Urban") +
#  geom_ribbon(predicted_df, mapping=aes(x=x, ymax=max.conf.high, ymin=min.conf.low, group=group), alpha=0.1)+
#  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'), legend.title=element_blank())
#results.plot
#
#ggplot()+
##  geom_point(dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.5, alpha=0.2)+
#  geom_line(predicted.full, mapping=aes(x=x, y=predicted, color=group))+
#  geom_ribbon(predicted_df, mapping=aes(x=x, ymax=max.conf.high, ymin=min.conf.low, group=group), alpha=0.1)
## this definitely does not work









################## 4. Trying spatial models in the spatialreg package ###################


# sample data
set.seed(10)
dat.samp <- dat[sample(nrow(dat), 5000), ]
# turn into sf object
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtered.tif")
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat"), crs=st_crs(GHSL)) 

mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat.samp) # here is the model that I am working with
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
                        precip + log(number_checklists) + elevation, data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
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
                                       precip + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
beep()
summary(dat.samp.slm)

# Spatial error model
dat.samp.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                         precip + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
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
                     precip + log(number_checklists), data = dat.samp2)

# turn into sf object
dat.samp2.sf <- st_as_sf(dat.samp2, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.samp2.nb <- dnearneigh(dat.samp2.sf, d1=0, d2=5)
dat.samp2.lw <- nb2listw(dat.samp2.nb, style = "W", zero.policy = TRUE)

dat.samp2.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                          precip + log(number_checklists), data = dat.samp2, listw = dat.samp2.lw, zero.policy = TRUE)

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
                                    precip + log(number_checklists), data = dat, listw = dat.lw, zero.policy = TRUE)
# can't run, sample size too large
dat.samp2$residuals.lm <- residuals(dat.samp2.lm)
moran.mc(dat.samp2$residuals.lm, dat.samp2.lw, nsim = 999, zero.policy = TRUE)














#########################################################

##### Thinning of proportional loss in richness
dat <- read.csv("modeling_data.csv")

# Run model
full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                precip + log(number_checklists) + elevation, dat)





predicted.Nhemisphere <- predictions(
  full.model,
  type = "response", transform=square,
  # by = "urban2",
  newdata = datagridcf(urban2=c("Natural", "Urban", "Suburban"))) %>% 
  select(estimate, abslat, urban2, hemisphere) %>% filter(hemisphere=="northern") %>% pivot_wider(names_from="urban2", values_from=c("estimate"), values_fn=mean)

# S hemisphere
predicted.Shemisphere <- predictions(
  full.model,
  type = "response", transform=square,
  # by = "urban2",
  newdata = datagridcf(urban2=c("Natural", "Urban", "Suburban"))) %>% 
  select(estimate, abslat, urban2, hemisphere) %>% filter(hemisphere=="southern") %>% pivot_wider(names_from="urban2", values_from="estimate", values_fn=mean)

predicted.fulldata <- rbind(predicted.Nhemisphere, predicted.Shemisphere)

# merge with observed values
predicted.fulldata2 <- inner_join(total.dat[, c(2, 3, 21, 26, 27, 28, 31, 32, 34)], predicted.fulldata, by=c("abslat", "hemisphere"))

# line for urban sites
resids <- predicted.fulldata2 %>% group_by(urban2) %>% mutate(resids = total_SR/Natural)

# plot
ggplot(resids, aes(x=abslat, y=resids, group=urban2, color=urban2))+
  geom_point(size=0.2, alpha=0.2)+
  facet_wrap(~hemisphere)+
  geom_smooth(method="lm")

# now I want to model the residuals as a linear model

hist(sqrt(resids$resids)) # they have a square root distribution, which makes sense because the data they are pulled from has a square root dist

mod.resids <- lm(sqrt(resids)~abslat * urban2 * hemisphere, resids)
plot(mod.resids)
# looks pretty good actually

# Now plot model predictions
predicted.resids <- avg_predictions(mod.resids, by=c("abslat", "urban2"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))

ggplot()+
  geom_point(resids, mapping=aes(x=abslat, y=resids, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.resids, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.resids, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Proportion richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15))
# ok that looks pretty good


# but then this is probably autocorrelated as well so I need to do the thinning
mod.resids.gls <- gls(sqrt(resids)~abslat * urban2 * hemisphere, resids)
library(sf)
resids.sf <- st_as_sf(resids, coords=c("long", "lat")) 

plot(gstat::variogram(residuals(mod.resids.gls, "normalized") ~
                        1, data = resids.sf, cutoff = 200))
# doesn't look like there is much autocorrelation
library(spdep)

# there is autocorrelation on sample so try doing the thinning
library(terra)
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL))

resids.samp <- resids[sample(nrow(resids), 10000), ]

vect <- st_as_sf(resids.samp, crs=st_crs(4326), coords=c("long","lat"))
vect2 <- st_transform(vect, crs=crs(GHSL))

xy=st_coordinates(vect2)
# get cell number that each point is in
resids.samp$cell.subsample<-cellFromXY(sample.grid, xy)

# randomly sample one point within each cell
resids.thinned <- resids.samp %>% group_by(cell.subsample) %>% sample_n(1) 



lm.resids.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                       precip + log(number_checklists) + elevation, resids.thinned)
resids.thinned.sf <- st_as_sf(resids.thinned, coords=c("long", "lat")) 
resids.thinned.nb <- dnearneigh(resids.thinned.sf, d1=0, d2=200) # calculate distances
resids.thinned.lw <- nb2listw(resids.thinned.nb, style = "W", zero.policy = TRUE)
## Run Moran's I
moran.samp <- lm.morantest(lm.resids.thinned, resids.thinned.lw, zero.policy = T)
moran.samp
library(beep)
beep()
















