## In this script, I will start running some models of how the latitudinal diversity gradient changes with urbanization
library(terra)
library(sf)
library(tidyverse)
library(nlme)
library(ncf)
library(automap)
library(spdep)
library(gstat)

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





######################################
# Start running models

dat <- read.csv("modeling_data.csv")
summary(dat)
hist(dat$total_SR) # a bit skewed, also count data
hist(log(dat$total_SR))
hist(sqrt(dat$total_SR)) # this looks pretty good
hist(dat$number_checklists) # this is super log normal, used the log in the response variable
# Try a simple linear model with absolute latitude
mod1 <- lm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + BIOME + log(number_checklists), dat)
mod1.trans <- lm(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

summary(mod1)
summary(mod1.trans)
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
 #                          z=dat$gls1Resids, resamp=100, quiet=TRUE) # this takes up too much memory in R because there are so many data points

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
#w <-as.matrix(1/dist(cbind(dat$long, dat$lat))) # this is too much for my computer to handle, do on lab computer
#wlist<-mat2listw(w)
#moran.test(dat$gls1Resids, wlist)

# Fit all possible autocorrelation structures and see which is the best fit
csSpher <- corSpher(form=~lat+long,nugget=TRUE) # sperical
csExp <- corExp(form=~lat+long,nugget=TRUE) # exponential
csGaus <- corGaus(form=~lat+long,nugget=TRUE) # gaussian
csLin <- corLin(form=~lat+long,nugget=TRUE) # linear
csRatio <- corRatio(form=~lat+long,nugget=TRUE) # ratio
# update models
glsSpher <- update(gls1, correlation=corSpher(form=~lat+long,nugget=TRUE)) # spherical
glsExp <- update(gls1, correlation=csLin)

############### Subsample to test things
# since some things are not working, going to see if it is a problem with the sample size (and my computer memory)
# going to subsample within my data
dat.samp <- dat[sample(nrow(dat), 5000), ]
# try to run a correlog


gls1.samp <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + BIOME + log(number_checklists), dat.samp)

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
# it is signficant


lsSpher <- update(gls1.samp, correlation=corSpher(form=~lat+long,nugget=TRUE)) # spherical










