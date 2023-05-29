########### This is the script for modelling the seasonal differences in bird diversity ###################
library(tidyverse)
library(terra)
library(sp)
library(nlme)
library(ncf)

# Load data
dat <- read.csv("season_model_data.csv")

# model
mod1 <- lm(sqrt(total_SR) ~ abs(lat) * urban * season + hemisphere + CONTINENT +
                   abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

summary(mod1)
anova(mod1)
# triple interaction is significant

# look at model residuals
plot(mod1)
# looking pretty good

# do another one with absolute lat for plotting
dat$abslat <- abs(dat$lat)
mod1.abslat <- lm(sqrt(total_SR) ~ abslat * urban * season + hemisphere + CONTINENT +
             abs(lat):CONTINENT + BIOME + log(number_checklists), dat)

# run a gls
gls1 <- gls(sqrt(total_SR) ~ abs(lat) * urban * season + hemisphere + CONTINENT + BIOME + log(number_checklists), dat)

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
library(interactions)
interact_plot(mod1, lat, season)
# steeper gradient in winter as predicted

# try to make a correlogram
dat$gls1.resids <- residuals(gls1)
hist(dat$gls1.resids) # looking pretty good
residsI <- spline.correlog(x=dat$long, y=dat$lat, z=dat$gls1.resids, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20)) 




