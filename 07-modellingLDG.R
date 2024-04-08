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
library(marginaleffects)
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

dat %>% group_by(precip) %>% summarise(n=n()) # look at how many observations per precip

dat %>% group_by(urban2) %>% count()

ggplot(dat, aes(x=abslat, y=log(number_checklists), color=urban2))+
  geom_point(alpha=0.2)+
  geom_smooth()










###################################
# START MODELLING


dat %>% group_by(CONTINENT) %>% summarise(n=n()) # most are in quadrant 2 which includes north america
# Try a simple linear model with absolute latitude

mod1 <- lm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
              abs(lat):CONTINENT + precip + log(number_checklists), dat)



# try a bunch of different models
mod1.trans <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat) # latitude and hemisphere interaction

mod1.hemsiphere <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + 
                   precip + log(number_checklists) + elevation, dat) # no interaction with hemisphere

mod1.trans.cont <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + abslat:CONTINENT + 
                        precip + log(number_checklists) + elevation, dat) # continent as a fixed effect


mod1.trans.cont.intrxn <- lm(sqrt(total_SR) ~ abslat * urban2 * CONTINENT + 
                             precip + log(number_checklists) + elevation, dat) # triple interaction between continent, latitude, and urbanization


AIC(mod1.trans, mod1.hemsiphere, mod1.trans.cont, mod1.trans.cont.intrxn)
# model with continent interaction is best, but I don't think I want to use this one
# I think I will do the hemisphere one with the interaction

#mod1.quadrant <- lm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
 #                     precip + log(number_checklists) + elevation, dat) # model with quadrant instead

###### Model with and without interaction between everything and number of checklists
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat)
plot_slopes(mod1, variables="abslat", condition=c("urban2")) # all very much different
summary(mod1)

# compute partial R2 and partial cohen's F2
library(sensemakr)
partial_r2 <- as.data.frame(partial_r2(mod1))
partial_f2 <- as.data.frame(partial_f2(mod1))





mod2 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere * log(number_checklists) +
             precip + elevation, dat)
summary(mod2)
anova(mod2)
# quadruple interaction is not significant, can probably remove
AIC(mod1, mod2) # mod2 is better
anova(mod2, mod1) # mod2 is better


# The quadruple interaction is a bit complicated, let's try with simpler
mod3 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + log(number_checklists) + 
             hemisphere:log(number_checklists) + urban2:log(number_checklists) + abslat:log(number_checklists) +
             precip + elevation, dat)

plot_slopes(mod3, variables="abslat", condition=c("urban2")) # all very much different
plot_slopes(mod2, variables="abslat", condition=c("urban2")) # all very much different



plot_predictions(mod1, by=c("abslat", "urban2"), transform=square, 
                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))
plot_predictions(mod2, by=c("abslat", "urban2"), transform=square, 
                 newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))
plot_predictions(mod3, by=c("abslat", "urban2"), transform=square, 
                 newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))


######## Look at results
means.df <- as.data.frame(emmeans(mod1.trans, specs="urban2"))
means.df$emmean.sq <- means.df$emmean^2
135-112
111-93

emtrends(mod1.trans, pairwise ~ urban2, var="abslat")
summary(mod1.trans)

emtrends(mod1.trans, pairwise ~ urban2, var="abslat", by="hemisphere")
# all significantly negative

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
predicted3 <- ggpredict(mod1.trans, terms = c("abslat", "urban2", "hemisphere")) 
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

effect_plot(mod1.trans, pred=precip, interval=TRUE) # generally a positive relationship between precipitation and richness

effect_plot(mod1.trans, pred=CONTINENT, interval=TRUE)


library(interactions)


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
                 lat:CONTINENT + precip + log(number_checklists), dat)
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
                   abs(lat):CONTINENT + precip + log(number_checklists), dat) # looks good
AIC(mod1.trans)
AIC(mod.3cat) # less good of a model but not too far off
mod1.abslat <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + CONTINENT +
                    abs(lat):CONTINENT + precip + log(number_checklists), dat)

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
              abs(lat):CONTINENT + precip + log(number_checklists), dat.97)
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
                   abs(lat):CONTINENT + precip + log(number_checklists), dat.samp2)

glsExp2 <- update(gls2.samp, correlation=csExp)



############ Trying a poisson model
mod.poisson <- glm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
      abs(lat):CONTINENT + precip + log(number_checklists), data=dat, family=poisson)
plot(mod.poisson)
summary(mod.poisson) # it is overdispersed
# everything very significant
# plot
ggplot(dat, aes(x=abs(lat), y=total_SR)) + geom_point() +
  stat_smooth(method="glm", method.args=list(family="poisson"), se=FALSE)+facet_wrap(~urban)

# try poly model with latitude
mod.poisson.poly <- glm(total_SR ~ lat + urban + lat:urban + I(lat^2) + lat:urban +  I(lat^2):urban 
                   + hemisphere + CONTINENT + precip + log(number_checklists), data=dat, family=poisson)

AIC(mod.poisson, mod.poisson.poly)
# the poly is worse




############################




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
                   precip + log(number_checklists) + elevation, data = samps) # with hemisphere

summary(dat.binned.mod1)
# everything still significant

dat.binned.mod2 <- lm(sqrt(total_SR) ~ abslat * urban2 + CONTINENT + 
                        precip + log(number_checklists) + elevation, data = samps) # with continent

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
                   precip + log(number_checklists) + elevation, data = samps, listw = samps.lw, zero.policy = TRUE)
summary(impacts(samps.slx, listw=samps.lw), zstats=TRUE)

# try spatial lag model
samps.slm <- lagsarlm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + abslat:hemisphere + 
                     precip + log(number_checklists), data = samps, listw = samps.lw, zero.policy = TRUE)
summary(samps.slm) # still spatially autocorrelated

# test for autocorrelation
samps$residuals.slx <- residuals(samps.slx)
moran.mc(samps$residuals.slx, samps.lw, nsim = 999, zero.policy = TRUE)










####################################
# TRY RUNNING SOME NON-LINEAR MODELS ON THE FULL DATASET TO SEE IF THERE IS A NON-LINEAR TREND
library(mgcv)
# model with regular latitude
mod.gam1 <- gam(total_SR ~ s(lat, by=c(urban2)) + hemisphere + precip + log(number_checklists) + elevation, data = dat)

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
                  precip + log(number_checklists) + elevation, data = dat) # couldn't do hemisphere intrxn for some reason

plot(ggeffects::ggpredict(mod.gam2, terms=c("abslat", "urban2"), facets = TRUE))
plot(ggeffects::ggpredict(mod.gam2, terms=c("elevation"), facets = TRUE))
# SR generally decreases by elevation
# this is the relationship with just absolute latitude overall
# there is not a hump in the middle even though there are more points there

# Try model without smoothing parameter and see how it compares
mod.gam3 <- gam(total_SR ~ abslat * urban2 + CONTINENT +
                  precip + log(number_checklists) + elevation, data = dat)

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
                  precip + log(number_checklists) + elevation, data=dat.n)
plot(ggeffects::ggpredict(mod.gam4, terms=c("lat", "urban2"), facets = TRUE)) # linear terms
# not a hump in the middle, that's good

dat.s <- dat %>% filter(hemisphere=="southern") 

mod.gam5 <- gam(total_SR ~ s(lat, by=urban2) +
                  precip + log(number_checklists) + elevation, data=dat.s)
plot(ggeffects::ggpredict(mod.gam5, terms=c("lat", "urban2"), facets = TRUE))
# who southern pattern is crazy but still not hump, that is good














###################################################
### Calculate prediction interval for my model
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2, dat) 

# creating new data
abslat <- sample(1:90, 60000, replace=TRUE)
urban2 <- as.factor(sample(c('Natural', 'Suburban', 'Urban'), 60000, replace=TRUE))
hemisphere <- sample(c('northern', 'southern'), 60000, replace=TRUE)
precip <- as.factor(sample(1:14, 60000, replace=TRUE))
number_checklists <- sample(4:11, 60000, replace=TRUE)
elevation <- sample(-900:4900, 60000, replace=TRUE)
newdata <- data.frame(abslat, urban2, hemisphere, precip, number_checklists, elevation)
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




