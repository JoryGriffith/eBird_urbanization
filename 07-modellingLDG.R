## In this script, I will start running some models of how the latitudinal diversity gradient changes with urbanization
library(terra)
library(sf)
library(tidyverse)
library(nlme)
library(ncf)
library(automap)
library(spdep)
library(gstat)
library(beepr)

#world <- ne_countries(scale = "medium", type="map_units", returnclass = "sf")


# load thresholded summary data
dat <- read.csv("5yr_summary/summary_thresholded.csv")
dat <- rename(dat, lat = y, long = x)

# First I need to get the other data that I want to include in the model

######## Assign hemisphere
dat <- dat %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))

######### Assign continent
#dat_sf <- st_as_sf(dat, coords=c('long', "lat"), crs=st_crs(world))

#joined <- st_join(dat_sf, world)


#dat_joined <- as.data.frame(joined[,c(1:22,41)] %>% mutate(long = sf::st_coordinates(.)[,1],
                                                  #         lat = sf::st_coordinates(.)[,2]))# just keep country name
# extract continent using country name
#dat_joined$continent <- countrycode(sourcevar = dat_joined[,"name_long"],
 #                                    origin = "country.name",
#                                     destination = "continent")


##################
# Trying new way to do continent
continents <- st_read("/Volumes/Expansion/eBird/continent-poly/Continents.shp")
#plot(continents)

dat_sf <- st_as_sf(dat, coords=c('long', "lat"), crs=st_crs(continents))

dat_cont <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature) # joining by nearest feature


#########################
# Extract biome
# Classifying points into biomes using a terrestrial biomes shapefile
biomes <- st_read("/Volumes/Expansion/eBird/wwf_biomes/wwf_terr_ecos.shp")
class(biomes) # sf and data frame
# look at how many biomes there are
length(unique(biomes$REALM)) # 9 of these
length(unique(biomes$BIOME)) # 16 biomes
length(unique(biomes$ECO_NAME)) # 827 ecoregion names

# plot biome
#plot(biomes["BIOME"])

# want to extract biomes
dat_withbiome <- st_join(dat_cont, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)


# create seperate columns for lat long again
datFINAL <- as.data.frame(dat_withbiome[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
                                             lat = sf::st_coordinates(.)[,2]))

summary(datFINAL)
# save as csv
write_csv(datFINAL, "modeling_data.csv")




#####################################
######################################
# Start running models

dat <- read.csv("modeling_data.csv")
dat %>% group_by(urban) %>% summarise(n=n())
summary(dat)
hist(dat$total_SR) # a bit skewed, also count data
hist(log(dat$total_SR))
hist(sqrt(dat$total_SR)) # this looks pretty good
hist(dat$number_checklists) # this is super log normal, used the log in the response variable

# turn biome and urban into a factor
class(dat$BIOME)
dat$BIOME <- as.factor(dat$BIOME)
dat$urban <- as.factor(dat$urban)
# Try a simple linear model with absolute latitude

mod1 <- lm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

mod1.trans <- lm(sqrt(total_SR) ~ abs(lat) * urban * hemisphere + CONTINENT +
                   abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

dat$abslat <- abs(dat$lat)
mod1.abslat <- lm(sqrt(total_SR) ~ abslat * urban * hemisphere + CONTINENT +
                   abs(lat):CONTINENT + BIOME + log(number_checklists), dat)
summary(mod1)
summary(mod1.trans)
anova(mod1.trans)
AIC(mod1)
AIC(mod1.trans) # the one with the sqrt tranformation is much lower
# these are looking pretty good
plot(mod1.trans)
plot(mod1) # sqrt transformation is a better fit and more normal
hist(residuals(mod1))
hist(residuals(mod1.trans)) # they are both pretty normally distributed but the transformed one is better


# rerun using gls (same as linear model when no spatial autocorrelation included)
gls1 <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

anova(gls1) # everything very significant

# check for spatial autocorrelation
dat$gls1Resids <- residuals(gls1)
hist(dat$gls1Resids) # looking very normally distributed, that's good

#residsI <- spline.correlog(x=dat$long, y=dat$lat,
                     #      z=dat$gls1Resids, resamp=100, quiet=TRUE) # this takes up too much memory in R because there are so many data points

#plot(residsI,xlim=c(0,50))



#### Variogram
#plot(nlme:::Variogram(gls1, form = ~lat +
#                        long, resType = "normalized"))
# this uses too much memory

# try a different function
datSPDF <- dat
coordinates(datSPDF) <- c("long","lat")
plot(gstat::variogram(residuals(gls1, "normalized") ~
                 1, data = datSPDF, cutoff = 100))

# try with directional
plot(gstat::variogram(residuals(gls1, "normalized") ~
                        1, data = datSPDF, cutoff = 100, alpha = c(0, 45, 90, 135)))
# it likely looks like this at 0 because it becomes more similar as you move toward the equator again 

#### Moran's I
# make distance matrix
#w <-as.matrix(1/dist(cbind(dat$long, dat$lat))) # this is too much for my computer or the lab computer to handle
#wlist<-mat2listw(w)
#moran.test(dat$gls1Resids, wlist)

# Fit all possible autocorrelation structures and see which is the best fit
csSpher <- corSpher(form=~lat+long,nugget=TRUE) # sperical
csExp <- corExp(form=~lat+long,nugget=TRUE) # exponential
csGaus <- corGaus(form=~lat+long,nugget=TRUE) # gaussian
csLin <- corLin(form=~lat+long,nugget=TRUE) # linear
csRatio <- corRatio(form=~lat+long,nugget=TRUE) # ratio
# update models
glsGaus <- update(gls1, correlation=csLin) # gaussian
#glsExp <- update(gls1, correlation=csLin)



#####################
# Plotting model results
library(jtools)
effect_plot(mod1.trans, pred=urban, interval=TRUE) # negative relationship

effect_plot(mod1.trans, pred=lat, interval=TRUE) # peaks at intermediate latitudes, looks pretty good

effect_plot(mod1.trans, pred=hemisphere, interval=TRUE) # no difference

effect_plot(mod1.trans, pred=BIOME, interval=TRUE)

effect_plot(mod1.trans, pred=CONTINENT, interval=TRUE)


library(interactions)

interact_plot(mod1.trans, lat, hemisphere)
interact_plot(mod1.trans, lat, urban)
# woo this looks great!
# plot with absolute latitude
interact_plot(mod1.abslat, abslat, urban, plot.points=TRUE)
interact_plot(mod1.abslat, abslat, hemisphere)
#########################
# Going to try to do a quadratic model to see if it is a better fit
dat$lat2 <- (dat$lat)**2

head(dat)
mod.quad <- lm(sqrt(total_SR) ~ lat + lat2 + lat:urban + lat2:urban + hemisphere + CONTINENT +
                 lat:CONTINENT + BIOME + log(number_checklists), dat)
plot(mod.quad)
summary(mod.quad)
effect_plot(mod.quad, pred=list(lat,lat2), interval=TRUE)

AIC(mod.quad) # for some reason not showing up as quadratic
AIC(mod1.trans)
# the non-quadratic model is better
interact_plot(mod.quad, lat, urban)



# try model with only 3 categories
dat <- dat %>% mutate(urban2=ifelse(urban==11, 1, ifelse(urban==30, 3, 2)))
dat %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat$urban2 <- as.factor(dat$urban2)
mod.3cat <- lm(sqrt(total_SR) ~ abs(lat) * urban2 + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + BIOME + log(number_checklists), dat) # looks good
AIC(mod1.trans)
AIC(mod.3cat) # less good of a model but not too far off
mod1.abslat <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + CONTINENT +
                    abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

interact_plot(mod.3cat, lat, urban2, interval=TRUE)
interact_plot(mod1.abslat, abslat, urban2, interval=TRUE) 
# looking good!



############### Subsample to test things
# since some things are not working, going to see if it is a problem with the sample size (and my computer memory)
# going to subsample within my data
dat.samp <- dat[sample(nrow(dat), 5000), ]
# try to run a correlog

gls1.samp <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + as.factor(BIOME) + log(number_checklists), dat.samp)


# Make a correlogram
dat.samp$gls1.samp <- residuals(gls1.samp)
hist(dat.samp$gls1.samp)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$gls1.samp, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20)) # there is some autocorrelation at small distances


# calculate Moran's I
w <-as.matrix(1/dist(cbind(dat.samp$long, dat.samp$lat))) # make inverse distance matrix - weights things that are close together higher
wlist<-mat2listw(w) # assign weights based on inverse distances (converts square spatial weights matrix to a weights list object)
# this quickly becomes very large because it is pairwise distances
moran.test(dat.samp$gls1.samp, wlist)
# it is significant

# update with gaussian
glsGaus <- update(gls1.samp, correlation=csGaus) # gaussian
dat.samp$glsGaus <- residuals(glsGaus)
hist(dat.samp$glsGaus)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$glsGaus, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
beep()
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
plot(residsI,xlim=c(0,20))

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

# plot variograms
plot(nlme::Variogram(gls1.samp, form =~lat + long, resType="normalized")) # plot original model
plot(nlme::Variogram(glsSpher, form =~lat + long, resType="normalized")) # plot model with spherical autocorrelation structure
plot(nlme::Variogram(glsGaus, form =~lat + long, resType="normalized")) # gaussian
plot(nlme::Variogram(glsExp, form =~lat + long, resType="normalized")) # exponential
plot(nlme::Variogram(glsLin, form =~lat + long, resType="normalized")) # linear
plot(nlme::Variogram(glsSpher, form =~lat + long, resType="normalized")) # gaussian
# these all look very similar

# test moran's I
w <-as.matrix(1/dist(cbind(dat.samp$long, dat.samp$lat))) # make inverse distance matrix - weights things that are close together higher
wlist<-mat2listw(w) # assign weights based on inverse distances (converts square spatial weights matrix to a weights list object)
# this quickly becomes very large because it is pairwise distances
moran.test(dat.samp$glsSpher, wlist)
# all of the models with the spatial autocorrelation structures still have significant autocorrelation
# this is a problem


##########################################
# Try modelling with the other thresholds
# The 97% one is 139 checklists
# The 98% one is 203 checklists

dat.97 <- dat %>% filter(number_checklists >= 139) # 47057 observations

# Model 
gls2 <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + BIOME + log(number_checklists), dat.97)
anova(gls2) # all very significant

# try to add a correlation structure
glsExp2 <- update(gls2, correlation=csExp)
# still too large

## Testing how many datapoints can be in the model and I can still be able to add a correlation structure
dat.samp2 <- dat[sample(nrow(dat), 40000), ]
# try to run a correlog
# seems to work with 40000 but not 50000
# with the highest threshold we have 37000 so maybe I should use that

gls2.samp <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + as.factor(BIOME) + log(number_checklists), dat.samp2)

glsExp2 <- update(gls2.samp, correlation=csExp)







