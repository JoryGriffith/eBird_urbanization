## In this script, I will start running some models of how the latitudinal diversity gradient changes with urbanization

# Load packages
library(terra)
library(sf)
library(tidyverse)
library(nlme)
library(ncf)
library(automap)
library(spdep)
library(gstat)
library(beepr)
library(jtools)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggeffects)
library(foreach)
library(doParallel)
library(spatialreg)
library(elevatr)


######################################
# Start running models

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
######

# Try a simple linear model with absolute latitude

mod1 <- lm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + BIOME + log(number_checklists), dat)


# try a bumch of different models
mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists), dat) # latitude and hemisphere interaction

mod1.hemisphere.intrxn <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   BIOME + log(number_checklists), dat) # urban, latitude, and hemisphere triple interaction

mod1.trans.cont <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + abslat:CONTINENT + 
                        BIOME + log(number_checklists), dat) # continent and latitude interaction

mod1.trans.wele <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists) + elevation, dat) # hemisphere and latitude interaction with elevation

mod1.trans.cont.wele <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + abslat:CONTINENT + 
                        BIOME + log(number_checklists) + elevation, dat) # continent and latitude interaction with elevation

mod1.trans.cont.intrxn <- lm(sqrt(total_SR) ~ abslat * urban2 * CONTINENT + 
                             BIOME + log(number_checklists) + elevation, dat) # triple interaction between continent, latitude, and urbanization

summary(mod1.trans.wele)
AIC(mod1.trans, mod1.hemisphere.intrxn, mod1.trans.cont, mod1.trans.wele, mod1.trans.cont.wele, mod1.trans.cont.intrxn) 
# last one with continent urbanization latitude interaction is best

# Plot model results for talk
predicted <- ggpredict(mod1.trans.cont.intrxn, terms = c("abslat", "urban2")) 
# looks the same whether sqrt included in model or not


results.plot <-
  plot(predicted, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

#ggsave(results.plot, file="results.plot.png", height=5, width=9)


# Compare slopes
lstrends(mod1.trans, pairwise ~ urban2, var="abslat")


# Compare between continents
predicted2 <- ggpredict(mod1.trans.cont.intrxn, terms = c("abslat", "urban2", "CONTINENT")) 
# looks the same whether sqrt included in model or not

results.plot2 <-
  plot(predicted2, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
results.plot2
ggsave(results.plot2, file="results.bycontient.png")

# Compare between hemispheres
predicted3 <- ggpredict(mod1.hemisphere.intrxn, terms = c("abslat", "urban2", "hemisphere")) 
# looks the same whether sqrt included in model or not

results.plot3 <-
  plot(predicted3, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
results.plot3 # gradient steeper in the southern hemisphere
ggsave(results.plot3, file="results.hemisphere.png")


# Plotting model results using jtools and interactions packages
library(jtools)
effect_plot(mod1.trans, pred=urban, interval=TRUE) # negative relationship

effect_plot(mod1.trans, pred=lat, interval=TRUE) # peaks at intermediate latitudes, looks pretty good

effect_plot(mod1.trans, pred=hemisphere, interval=TRUE) # no difference

effect_plot(mod1.trans, pred=BIOME, interval=TRUE)

effect_plot(mod1.trans, pred=CONTINENT, interval=TRUE)


library(interactions)

interact_plot(mod1.trans, lat, hemisphere)

interact_plot(mod1.trans, abslat, urban, interval=TRUE)

# woo this looks great!
# plot with absolute latitude
interact_plot(mod1.trans, abslat, urban)
interact_plot(mod1.abslat, abslat, hemisphere)

###########################
# LOOK AT SPATIAL AUTOCORRELATION

# rerun using gls (same as linear model when no spatial autocorrelation included)
gls1 <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + abs(lat):hemisphere + 
              BIOME + log(number_checklists), dat)

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
# This takes forever to run



#####################
#########################
# Going to try to do a quadratic model to see if it is a better fit

head(dat)
mod.quad <- lm(sqrt(total_SR) ~ lat + urban + lat:urban + I(lat^2) + lat:urban +  I(lat^2):urban 
               + hemisphere + CONTINENT +
                 lat:CONTINENT + BIOME + log(number_checklists), dat)
plot(mod.quad)
summary(mod.quad)
effect_plot(mod.quad, pred=lat, interval=TRUE)
interact_plot(mod.quad, lat, urban)

AIC(mod.quad) # for some reason not showing up as quadratic
AIC(mod1.trans)

# the non-quadratic model is better
interact_plot(mod.quad, lat, urban)
# this looks weird, not sure what is happening


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

# need to look into the rank deficient fit error - says it may because there are colinear covariates or having more parameters than available variables





##########################
### SUBSAMPLING

# since some things are not working, going to see if it is a problem with the sample size (and my computer memory)
# going to subsample within my data
dat.samp <- dat[sample(nrow(dat), 1000), ]
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



############ Trying a poisson model
mod.poisson <- glm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
      abs(lat):CONTINENT + as.factor(BIOME) + log(number_checklists), data=dat, family=poisson)
plot(mod.poisson)
summary(mod.poisson) # it is overdispersed
# everything very significant
# plot
ggplot(dat, aes(x=abs(lat), y=total_SR)) + geom_point() +
  stat_smooth(method="glm", method.args=list(family="poisson"), se=FALSE)+facet_wrap(~urban)

# try poly model with latitude
mod.poisson.poly <- glm(total_SR ~ lat + urban + lat:urban + I(lat^2) + lat:urban +  I(lat^2):urban 
                   + hemisphere + CONTINENT + as.factor(BIOME) + log(number_checklists), data=dat, family=poisson)

AIC(mod.poisson, mod.poisson.poly)
# the poly is worse




############################
# TRYING THINNING

# look at the spread of the data between continents
dat %>% group_by(CONTINENT) %>% summarise(n=n())
# there are 55,000 in North America compared to much less everywhere else
# I could try to subsample in just north america and see if that reduces the spatial autocorrelation

# right now there are 70,749 points, lets see if I thinned that to 50,000
library(spThin)
# I will try to run spThin to a distance where there is no spatial autocorrelation, based on the correlogram of the sample data
# it looks like the autocorrelation ends at about 5 km so I am going to try thinning by that and see what happens
dat.samp <- dat[sample(nrow(dat), 1000), ]

dat.thinned <- thin.algorithm(dat.samp[,c(25:26)], thin.par=10, 
               rep=10)

da.thinned.int <- dat.thinned[[1]] %>% rename(long=Longitude, lat=Latitude)

dat.thinned.new <- inner_join(dat.samp, da.thinned.int, by=c("long", "lat"))

# run model with thinned data and check for autocorrelation
gls1.thinned <- gls(sqrt(total_SR) ~ abs(lat) * urban + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + as.factor(BIOME) + log(number_checklists), dat.thinned.new)
dat.thinned.new$gls1.thinned <- residuals(gls1.thinned)
hist(dat.thinned.new$gls1.thinned)
residsI <- spline.correlog(x=dat.thinned.new$long, y=dat.thinned.new$lat, z=dat.thinned.new$gls1.thinned, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
beep()

# test Moran's I
w <-as.matrix(1/dist(cbind(dat.thinned.new$long, dat.thinned.new$lat))) # make inverse distance matrix - weights things that are close together higher
?dist
wlist<-mat2listw(w) # assign weights based on inverse distances (converts square spatial weights matrix to a weights list object)
# this quickly becomes very large because it is pairwise distances
moran.test(dat.thinned.new$gls1.thinned, wlist)
beep()
# still spatially autocorrelated


# what if I add a correlation structure
glsSpher <- update(gls1.thinned, correlation=csSpher)
dat.thinned.new$glsSpher <- residuals(glsSpher)
hist(dat.thinned.new$glsSpher)
residsI <- spline.correlog(x=dat.thinned.new$long, y=dat.thinned.new$lat, z=dat.thinned.new$glsSpher, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
beep()

w <-as.matrix(1/dist(cbind(dat.thinned.new$long, dat.thinned.new$lat))) # make inverse distance matrix - weights things that are close together higher
wlist<-mat2listw(w) # assign weights based on inverse distances (converts square spatial weights matrix to a weights list object)
# this quickly becomes very large because it is pairwise distances
moran.test(dat.thinned.new$glsSpher, wlist)
beep()
# this is not working for some reason - all the residuals are the same

# try another correlation structure
glsGaus <- update(gls1.thinned, correlation=csGaus)
dat.thinned.new$glsGaus <- residuals(glsGaus)
hist(dat.thinned.new$glsGaus)
residsI <- spline.correlog(x=dat.thinned.new$long, y=dat.thinned.new$lat, z=dat.thinned.new$glsGaus, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
beep()

glsExp <- update(gls1.thinned, correlation=csExp)
dat.thinned.new$glsExp <- residuals(glsExp)
hist(dat.thinned.new$glsGaus)
residsI <- spline.correlog(x=dat.thinned.new$long, y=dat.thinned.new$lat, z=dat.thinned.new$glsExp, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
beep()



# look at variogram with and without correlation structure
plot(nlme::Variogram(gls1.thinned, form =~lat + long, resType="normalized"))
plot(nlme::Variogram(glsSpher, form =~lat + long, resType="normalized"))
AIC(gls1.thinned, glsSpher, glsExp)
# better with the correlation structure
beep()

anova(gls1.thinned)
anova(glsSpher)




###################################
#### Trying new ways to deal with spatial autocorrelation using the spatialreg package

# sample data
set.seed(10)
dat.samp <- dat[sample(nrow(dat), 10000), ]
# turn into sf object
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 

mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists), dat.samp) # here is the model that I am working with
summary(mod1.trans)
dat.samp$residuals <- residuals(mod1.trans)
dat.samp$fitted <- fitted(mod1.trans)
#saveRDS(mod1.trans, "50K_samp_lmmod.rds")

dat.samp.nb <- dnearneigh(dat.samp.sf, d1=0, d2=100) # calculate distances
# This function identifies neighbours of region points by euclidean distance
# it returns a list of integer vectors giving the region id numbers for neighbors satisfying the distance criteria
class(dat.samp.nb) # this is an nb object

?dnearneigh
dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)

# Moran's I test
lm.morantest(mod1.trans, dat.samp.lw, zero.policy = T) # very spatially autocorrelated
beep()
lm.LMtests(mod1.trans, dat.samp.lw, test="all", zero.policy = T) # test for spatial error - very significant
beep()
#lm.LMtests(mod1.trans, dat.samp.lw, test="LMlag", zero.policy = T) # test for spatial lag - also very significant
# LMerr and LMlag both significant, RLMlag significant, SARMA significant


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

for (d in seq(50, 2000, 50)) {
#foreach  (d = seq(50, 2000, 50),
 #                      .combine = 'c') %dopar% {
  dat.samp.nb <- dnearneigh(dat.samp.sf, d1 = 0, d2 = d)
  dat.samp.lw <- nb2listw(dat.samp.nb, style = "W", zero.policy = TRUE)
  moran <- moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE)
  moran_I <- c(moran_I, moran$statistic)
} # THIS TAKES A REALLY LONG TIME

moran_I <- data.frame(moran = moran_I, 
                      distance = seq(50, 2000, 50))

ggplot(moran_I, aes(x = distance, y = moran)) + 
  geom_point() +
  geom_line()
beep()



##### Try running models  with sample data
# Spatially lagged X model
dat.samp.slx <- lmSLX(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                                       BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
summary(dat.samp.slx) # R2 is 0.38 (higher than without the lag)
summary(mod1.trans)
summary(impacts(dat.samp.slx, listw=dat.samp.lw), zstats=TRUE) # indirect impacts are still pretty high
#saveRDS(dat.samp.slx, "50K_samp_slxmod.rds")

# test model
dat.samp$residuals.slx <- residuals(dat.samp.slx)
moran.mc(dat.samp$residuals.slx, dat.samp.lw, nsim = 999, zero.policy = TRUE)




# Spatial lag model
dat.samp.slm <- spatialreg::lagsarlm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                       BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
summary(dat.samp.slm)

# Spatial error model
dat.samp.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
           BIOME + log(number_checklists), data = dat.samp, listw = dat.samp.lw, zero.policy = TRUE)
summary(dat.samp.sem)



AIC(dat.samp.slx, dat.samp.sem, dat.samp.slm)
# the slx model is best

# run moran test
dat.samp$residuals.slx <- residuals(dat.samp.slx)
dat.samp$residuals.slm <- residuals(dat.samp.slm)
dat.samp$residuals.sem <- residuals(dat.samp.sem)

moran.mc(dat.samp$residuals, dat.samp.lw, nsim = 999, zero.policy = TRUE) # significant p = 0.001
moran.mc(dat.samp$residuals.slx, dat.samp.lw, nsim = 999, zero.policy = TRUE) # this is not significant!!! This might be the best model!
moran.mc(dat.samp$residuals.slm, dat.samp.lw, nsim = 999, zero.policy = TRUE) # significant but less p = 0.005
moran.mc(dat.samp$residuals.sem, dat.samp.lw, nsim = 999, zero.policy = TRUE) # significant p=0.001
# these all got rid of spatial autocorrelation!
beep()


# Load models with 50k samples
nonspatial.mod <- readRDS("50K_samp_lmmod.rds")
spatial.mod <- readRDS("50K_samp_slxmod.rds")

plot(ggeffects::ggpredict(spatial.mod, terms=c("urban2Suburban.n...17.623"), listw=dat.samp.lw, facets = TRUE))




# try to run on full model
dat.sf <- st_as_sf(dat, coords=c('long', "lat")) 

mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists), dat) # here is the model that I am working with
summary(mod1.trans)
dat$residuals <- residuals(mod1.trans)
dat$fitted <- fitted(mod1.trans)


dat.nb <- dnearneigh(dat.sf, d1=0, d2=100) # calculate distances (up to 200km)
dat.lw <- nb2listw(dat.nb, style = "W", zero.policy = TRUE)

dat.slx <- lmSLX(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                        BIOME + log(number_checklists), data = dat, listw = dat.lw, zero.policy = TRUE)
# DOESN'T WORK (TOO MANY SAMPLES)






############################################
# Bin the data by 5 degrees of latitude and sample within the bins
bin.size <- 5

chunks = split(dat, ceiling(dat$lat/bin.size)) # divide into bins
length(chunks)

datalist = vector("list", length = length(chunks))

set.seed(40)
for (i in 1:length(chunks)){
  if (nrow(chunks[[i]]) > 200) {
    chunki <- chunks[[i]]
  datalist[[i]] <- chunki[sample(nrow(chunki), 200), ] }
  else { (datalist[[i]] <- chunks[[i]])}
}

samps <- dplyr::bind_rows(datalist)

# Run model with more evenly sampled areas
dat.binned.mod1 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   BIOME + log(number_checklists) + elevation, data = samps) # with hemisphere

summary(dat.binned.mod1)
# everything still significant

dat.binned.mod2 <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + 
                        BIOME + log(number_checklists) + elevation, data = samps) # with continent

summary(dat.binned.mod2)

AIC(dat.binned.mod1, dat.binned.mod2) # better with continent

# Plot model results
predicted1 <- ggpredict(dat.binned.mod1, terms = c("abslat", "urban2")) 
# looks the same whether sqrt included in model or not

results.plot1 <-
  plot(predicted1, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
results.plot1

# Plot model results
predicted2 <- ggpredict(dat.binned.mod2, terms = c("abslat", "urban2")) 
# looks the same whether sqrt included in model or not

results.plot2 <-
  plot(predicted2, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
results.plot2 # this is a rank deficient fit because there are not enough points for an interaction for each continent
# looks different depending on whether I include continent as a fixed effect or an interaction
# But the basic pattern is the same

# test for spatial autocorrelation
samps.sf <- st_as_sf(samps, coords=c("long", "lat")) 
samps.nb <- dnearneigh(samps.sf, d1=0, d2=100) # calculate distances
samps.lw <- nb2listw(samps.nb, style = "W", zero.policy = TRUE) # turn into weighted list

# Moran's I test
lm.morantest(dat.binned.mod1, samps.lw, zero.policy = T)
lm.morantest(dat.binned.mod2, samps.lw, zero.policy = T)
# significant autocorrelation

# Try to run spatial model
samps.slx <- lmSLX(sqrt(total_SR) ~ abslat * urban2 * hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists), data = samps, listw = samps.lw, zero.policy = TRUE)
summary(impacts(samps.slx, listw=samps.lw), zstats=TRUE)

# try spatial lag model
samps.slm <- lagsarlm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + abslat:hemisphere + 
                     BIOME + log(number_checklists), data = samps, listw = samps.lw, zero.policy = TRUE)
summary(samps.slm)

plot(ggeffects::ggpredict(samps.slx, terms=c("abslat"), listw=samps.lw, facets = TRUE))

preds <- predict(samps.slx, listw=samps.lw, newdata=samps) # this does not work

# Trying to do this with OLS
x1 <- model.matrix(dat.binned.mod1)
lagx1 <- create_WX(x1,samps.lw,prefix="lagx")  # create lagged x variables
# Ok I am a bit confused what this is doing, says it is creating spatially lagged RHS variables
# But they are numerical when the variables are categorical, maybe it is based on weights? 
# Need to figure out exactly what this is doing
spat.data2 <- cbind(samps,lagx1) 
eg2c=lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + abslat:hemisphere + 
          BIOME + log(number_checklists) + lagx.SALESPC+ 
          lagx.COLLENRP +lagx.BKGRTOABC +lagx.BAPTISTSP +lagx.BKGRTOMIX +lagx.ENTRECP-1, data=spat.data2)



####################################
# TRY RUNNING SOME NON-LINEAR MODELS ON THE FULL DATASET TO SEE IF THERE IS A NON-LINEAR TREND
library(mgcv)

mod.gam1 <- gam(sqrt(total_SR) ~ s(abslat) * urban2 * hemisphere +
                  BIOME + log(number_checklists) + elevation, data = dat) # wrapping smoothing parameter around abslat
summary(mod.gam1)


predictions = predict(
  mod.gam1,
  newdata = dat,
  type = 'response',
  se = TRUE
)

df_preds = data.frame(dat, predictions) %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit)

ggplot(aes(x = abslat, y = total_SR, color=urban2), data = df_preds) +
  geom_smooth(method="gam")

plot(ggeffects::ggpredict(mod.gam1, terms = c("abslat", "urban2")), facets = TRUE) 
gratia::draw(mod.gam1)

# Try super simplified model
mod.gam2 <- gam(total_SR ~ s(abslat, by=urban2) + hemisphere +
                  BIOME + log(number_checklists) + s(elevation), data = dat) # couldn't do hemisphere intxn for some reason
summary(mod.gam2)
plot(ggeffects::ggpredict(mod.gam2, terms=c("abslat", "urban2"), facets = TRUE))
# SR generally decreases by elevation
# this is the relationship with just absolute latitude overall
# there is not a hump in the middle even though there are more points there


ggplot(dat, aes(y=total_SR, x=abslat, color=urban2)) +
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(), legend.spacing.y = unit(1, 'cm'))+
  geom_smooth(method="gam")
# there is not a peak
