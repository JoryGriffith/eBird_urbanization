#######################################
### This is the script for running the full model on compute canada
library(terra)
library(sf)
library(tidyverse)
library(spdep)
library(spatialreg)

# Load data for modelling
dat <- read.csv("modeling_data.csv")

# make data good for model
dat$BIOME <- as.factor(dat$BIOME)
dat$urban <- as.factor(dat$urban)

dat %>% group_by(urban) %>% summarise(n=n())
summary(dat)
hist(dat$total_SR, breaks=50)

dat$abslat <- abs(dat$lat)

dat %>% group_by(BIOME) %>% summarise(n=n())

# make another columbn with only 3 categories
# try model with only 3 categories
dat <- dat %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat$urban2 <- as.factor(dat$urban2)

dat$quadrant <- NA

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

dat <- dat %>% filter(!CONTINENT == "Antarctica")


# Run spatial model
GHSL <- rast("GHSL_filtered.tif")
dat.sf <- st_as_sf(dat, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.nb <- dnearneigh(dat.sf, d1=0, d2=5) # set distance to 5km
dat.lw <- nb2listw(dat.nb, style = "W", zero.policy = TRUE)

mod1.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                         BIOME + log(number_checklists), data = dat, listw = dat.lw, zero.policy = TRUE) # run spatial error model

saveRDS(mod1.sem, "spatialmod.rds")
saveRDS(dat.lw, "dat.lw.rds")


######################
## Run seasonal model
dat.seas <- read.csv("season_model_data.csv")
dat.seas$urban<-as.factor(dat.seas$urban)
dat.seas$BIOME <- as.factor(dat.seas$BIOME)

# make another column with just 3 categories
dat.seas <- dat.seas %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat.seas %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat.seas$urban2 <- as.factor(dat.seas$urban2)
dat.seas$abslat <- abs(dat.seas$lat) # absolute latitude

# Divide by quartiles
for (i in 1:nrow(dat.seas)){
  if (dat.seas$long[i] < 0 & dat.seas$hemisphere[i] == "northern") { # quadrant 1 is North America
    dat.seas$quadrant[i] <- 1
  }
  else if (dat.seas$long[i] > 0 & dat.seas$hemisphere[i] == "northern") { # quadrant 2 is europe and asia and N Africa
    dat.seas$quadrant[i] <- 2
  }
  else if (dat.seas$long[i] < 0 & dat.seas$hemisphere[i] == "southern") { # quadrant 3 is south america 
    dat.seas$quadrant[i] <- 3 
  }
  else {dat.seas$quadrant[i] <- 4} # quadrant 4 is oceania and southern africa
}
dat.seas$quadrant <- as.factor(dat.seas$quadrant)

dat.seas <- dat.seas %>% filter(!CONTINENT == "Antarctica") # filter out antarctica

## Spatial model
dat.seas.sf <- st_as_sf(dat.seas, coords=c("long", "lat"), crs=st_crs(GHSL)) 

dat.seas.nb <- dnearneigh(dat.seas.sf, d1=0, d2=5) # set distance to 5km
dat.seas.lw <- nb2listw(dat.seas.nb, style = "W", zero.policy = TRUE)

dat.seas.sem <- spatialreg::errorsarlm(sqrt(total_SR) ~ abslat * urban2 * quadrant + 
                                    BIOME + log(number_checklists), data = dat.seas, listw = dat.seas.lw, zero.policy = TRUE) # run spatial error model

saveRDS(dat.seas.sem, "spatialmod.seas.rds")
saveRDS(dat.seas.lw, "dat.seas.lw.rds")


















