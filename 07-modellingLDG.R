## In this script, I will start running some models of how the latitudinal diversity gradient changes with urbanization

# Load packages
library(terra)
library(sf)
library(tidyverse)
library(nlme)
library(ncf)
#library(automap)
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
library(emmeans)



############################################
# Load data
dat <- read.csv("modeling_data.csv")
dat %>% filter(CONTINENT=="Antarctica") # there is one point in antarctica, latitude is -54.05 so it just misses the cutoff
# this point is in the south sandwich islands

# look at everything above 65 degrees N
test <- dat %>% filter(lat >=65)
# there are 60 observations above 65 degrees N

hist(log(dat$total_SR))
hist(sqrt(dat$total_SR), breaks=50) # this looks pretty good
hist(dat$number_checklists) # this is super log normal, used the log in the response variable

dat %>% group_by(BIOME) %>% summarise(n=n()) # look at how many observations per biome





###################################
# START MODELLING


dat %>% group_by(CONTINENT) %>% summarise(n=n()) # most are in quadrant 2 which includes north america
# Try a simple linear model with absolute latitude

mod1 <- lm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + BIOME + log(number_checklists), dat)



# try a bunch of different models
mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   BIOME + log(number_checklists) + elevation, dat) # latitude and hemisphere interaction


mod1.trans.cont <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + abslat:CONTINENT + 
                        BIOME + log(number_checklists)+ elevation, dat) # continent and latitude interaction

mod1.trans.wele <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                   BIOME + log(number_checklists) + elevation, dat) # hemisphere and latitude interaction with elevation

mod1.trans.cont.wele <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + abslat:CONTINENT + 
                        BIOME + log(number_checklists) + elevation, dat) # continent and latitude interaction with elevation

mod1.trans.cont.intrxn <- lm(sqrt(total_SR) ~ abslat * urban2 * CONTINENT + 
                             BIOME + log(number_checklists) + elevation, dat) # triple interaction between continent, latitude, and urbanization






mod1.quadrant <- lm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                      BIOME + log(number_checklists) + elevation, dat) # model with quadrant instead


mod <-lm(sqrt(total_SR) ~ abslat * urban2, dat)

predicted <- ggpredict(mod1.trans, terms = c("abslat", "urban2")) 
# looks the same whether sqrt included in model or not

results.plot <-
  plot(predicted, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=15), legend.spacing.y = unit(1, 'cm'), legend.title=element_blank())
results.plot
# same results ! This is good
# save results plot
ggsave(results.plot, file="LDGMainResults.png", height=5, width=8)


######## Look at results
emmeans(mod1.trans, specs="urban2")
11.6^2 # 134.56
10.5^2 # 110
9.6^2 # 92
emtrends(mod1.trans, pairwise ~ urban2, var="abslat")
summary(mod1.trans)

emtrends(mod1.trans, pairwise ~ urban2, var="abslat", by="hemisphere")
# all significantly negative

summary(mod1.trans)
AIC(mod1.trans, mod1.hemisphere.intrxn, mod1.trans.cont, mod1.trans.wele, mod1.trans.cont.wele, mod1.trans.cont.intrxn, mod1.quadrant) 
# last one with continent urbanization latitude interaction is best
# The model with quadrant is better than the model with only hemisphere but not as good as the model with continent

# Plot model results for talk
predicted <- ggpredict(mod1.trans , terms = c("abslat", "urban2", "quadrant")) 
# looks the same whether sqrt included in model or not


results.plot <-
  plot(predicted, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
results.plot
ggsave(results.plot, file="results.plot.quadrant.png", height=5, width=9)


# Compare slopes
lstrends(mod1.trans, pairwise ~ urban2, var="abslat")


# Compare between continents
predicted2 <- ggpredict(mod1.trans.cont.intrxn, terms = c("abslat", "urban2")) 
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

###################################
#### USING SPATIALREG PACKAGE TO RUN MODELS

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
samps.nb <- dnearneigh(samps.sf, d1=0, d2=40) # calculate distances
samps.lw <- nb2listw(samps.nb, style = "W", zero.policy = TRUE) # turn into weighted list
# supplements a neighbors list with spatial weights for the chosen coding scheme
?nb2listw
# Moran's I test
lm.morantest(dat.binned.mod1, samps.lw, zero.policy = T)
lm.morantest(dat.binned.mod2, samps.lw, zero.policy = T) 

# significant autocorrelation

# Try to run spatial model
samps.slx <- lmSLX(sqrt(total_SR) ~ abslat * urban2 * hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists) + elevation, data = samps, listw = samps.lw, zero.policy = TRUE)
summary(impacts(samps.slx, listw=samps.lw), zstats=TRUE)

# try spatial lag model
samps.slm <- lagsarlm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + abslat:hemisphere + 
                     BIOME + log(number_checklists), data = samps, listw = samps.lw, zero.policy = TRUE)
summary(samps.slm) # still spatially autocorrelated

# test for autocorrelation
samps$residuals.slx <- residuals(samps.slx)
moran.mc(samps$residuals.slx, samps.lw, nsim = 999, zero.policy = TRUE)










####################################
# TRY RUNNING SOME NON-LINEAR MODELS ON THE FULL DATASET TO SEE IF THERE IS A NON-LINEAR TREND
library(mgcv)
# model with regular latitude
mod.gam1 <- gam(total_SR ~ s(lat, by=c(urban2)) + hemisphere + BIOME + log(number_checklists) + elevation, data = dat)

summary(mod.gam1)

plot(ggeffects::ggpredict(mod.gam1, terms = c("lat", "urban2")), facets = FALSE) 
#gratia::draw(mod.gam1)

#predictions = predict(
 # mod.gam1,
#  newdata = dat,
 # type = 'response',
#  se = TRUE
#)

#df_preds = data.frame(dat, predictions) %>%
 # mutate(lower = fit - 1.96 * se.fit,
  #       upper = fit + 1.96 * se.fit)

#ggplot(aes(x = abslat, y = total_SR, color=urban2), data = df_preds) +
 # geom_smooth(method="gam")

# model with absolute latitude
mod.gam2 <- gam(total_SR ~ s(abslat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data = dat) # couldn't do hemisphere intrxn for some reason

plot(ggeffects::ggpredict(mod.gam2, terms=c("abslat", "urban2"), facets = TRUE))
plot(ggeffects::ggpredict(mod.gam2, terms=c("elevation"), facets = TRUE))
# SR generally decreases by elevation
# this is the relationship with just absolute latitude overall
# there is not a hump in the middle even though there are more points there

# Try model without smoothing parameter and see how it compares
mod.gam3 <- gam(total_SR ~ abslat * urban2 + CONTINENT +
                  BIOME + log(number_checklists) + elevation, data = dat)

plot(ggeffects::ggpredict(mod.gam3, terms=c("abslat", "urban2"), facets = TRUE)) # linear terms

AIC(mod.gam2, mod.gam3) # the linear model is actually better

ggplot(dat, aes(y=total_SR, x=abslat, color=urban2)) +
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(), legend.spacing.y = unit(1, 'cm'))+
  geom_smooth(method="gam")
# there is not a peak

###### Look at the trends by hemisphere
dat.n <- dat %>% filter(hemisphere=="northern") 

mod.gam4 <- gam(total_SR ~ s(lat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data=dat.n)
plot(ggeffects::ggpredict(mod.gam4, terms=c("lat", "urban2"), facets = TRUE)) # linear terms
# not a hump in the middle, that's good

dat.s <- dat %>% filter(hemisphere=="southern") 

mod.gam5 <- gam(total_SR ~ s(lat, by=urban2) +
                  BIOME + log(number_checklists) + elevation, data=dat.s)
plot(ggeffects::ggpredict(mod.gam5, terms=c("lat", "urban2"), facets = TRUE))
# who southern pattern is crazy but still not hump, that is good














###################################################
### Calculate prediction interval for my model
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2, dat) 

# creating new data
abslat <- sample(1:90, 60000, replace=TRUE)
urban2 <- as.factor(sample(c('Natural', 'Suburban', 'Urban'), 60000, replace=TRUE))
hemisphere <- sample(c('northern', 'southern'), 60000, replace=TRUE)
BIOME <- as.factor(sample(1:14, 60000, replace=TRUE))
number_checklists <- sample(4:11, 60000, replace=TRUE)
elevation <- sample(-900:4900, 60000, replace=TRUE)
newdata <- data.frame(abslat, urban2, hemisphere, BIOME, number_checklists, elevation)
newdata2 <- data.frame(abslat, urban2)

prediction_interval <- predict(mod1, newdata=newdata2, interval="prediction")

confidence_interval <- predict(mod1, newdata=newdata2, interval="confidence") # also do confidence to compare

# Try with simpler model of just latitude and urbanization


# plot prediction intervals
prediction_df <- cbind(newdata2, prediction_interval)

ggplot(prediction_df, aes(x=abslat, y=fit))+
  geom_smooth(method="lm") +
  geom_line(data = prediction_df, aes(x = abslat, y = lwr), linetype = "dashed", color = "grey") +
  geom_line(data = prediction_df, aes(x = abslat, y = upr), linetype = "dashed", color = "grey") +
  facet_wrap(urban2)
# is it bad that they are so large and overlapping?

# plot confidence intervals
confidence_df <- cbind(newdata2, confidence_interval)

ggplot(confidence_df, aes(x=abslat, y=fit))+
  geom_smooth(method="lm") +
  geom_line(data = confidence_df, aes(x = abslat, y = lwr), linetype = "dashed", color = "grey") +
  geom_line(data = confidence_df, aes(x = abslat, y = upr), linetype = "dashed", color = "grey") +
  facet_wrap(urban2)
# these are sooooo much smaller

## I should try this with my covariates as random effects maybe? 
# I'm not sure what to do because these are pretty messed up




