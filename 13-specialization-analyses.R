########### Script for analyzing whether specialist species are being disproportionately lost at lower latitudes
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)
## First I want to thin the data so that for each cell, there is only one row for each species 

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4"
          , "r4c2", "r4c4") # take out r4c1 and r4c4 because they are not in the final modelling data


model.data <- read.csv("modeling_data.csv") 
# also going to load this because I want to filter for cells that are in the final dataset to save space
unique(model.data$square)
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

datalist.years <- list()
datalist.names <- list()


for (i in 16:16){ # come back to 16, 18
  for (j in 5:6) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", names[i], "_", years[j], "_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% model.data$cell) # filter for cells that are in the final dataset
    dat.filt$SCIENTIFIC.NAME<- as.character(dat.filt$SCIENTIFIC.NAME)
    dat.filt$OBSERVATION.DATE<- as.character(dat.filt$OBSERVATION.DATE)
    dat.filt$OBSERVER.ID<- as.character(dat.filt$OBSERVER.ID)
    dat.filt$SAMPLING.EVENT.IDENTIFIER <- as.character(dat.filt$SAMPLING.EVENT.IDENTIFIER)
    dat.filt$OBSERVATION.COUNT <- as.character(dat.filt$OBSERVATION.COUNT)
    dat.filt$GROUP.IDENTIFIER <- as.character(dat.filt$GROUP.IDENTIFIER)
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  dat_uniquesp$x <- xFromCell(GHSL, dat_uniquesp$cell) # extract the coordinates from the cells
  dat_uniquesp$y <- yFromCell(GHSL, dat_uniquesp$cell)
  dat_uniquesp$square=names[i]
  datalist.names[[i]] <- dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(dat_uniquesp)
}

global_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(global_unique_sp$SCIENTIFIC.NAME)) # 10,723 species
write.csv(global_unique_sp, paste("global_unique_species.csv"), row.names=FALSE)

# convert to lat long points to more easily bin by latitude
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
global_uniquesp <- read.csv("global_unique_species.csv")

dat_latlong <- st_as_sf(global_unique_sp, coords=c("x", "y"), crs=st_crs(GHSL))
dat_latlong <- st_transform(dat_latlong, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- as.data.frame(dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                     lat = sf::st_coordinates(.)[,2]))

# bind this with data in other crs
global_uniquesp <- cbind(global_unique_sp, latlong_df[,5:6])

#### Extract urban scores
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtMollweide.tif")
global_uniquesp$urban <- as.data.frame(terra::extract(GHSL, global_uniquesp[,c(3:4)]))$SMOD_global

### Merge with trait data
# Use taxise to make sure taxonomy lines up

diet<- read.csv("/Volumes/Backup/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv")

habitat <- read.csv("/Volumes/Backup/eBird/Traits/habitat_breadth.csv")

sp_diet <- merge(global_uniquesp, diet[, c(9,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific")
length(unique(sp_diet$SCIENTIFIC.NAME)) # 6,904 species

sp_habitat <- merge(global_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
length(unique(sp_habitat$SCIENTIFIC.NAME)) # 8,498 species

# Bin latitude by 10 degrees
sp_diet <- sp_diet %>% mutate(lat_bin = cut(lat, breaks=c(-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)))
# see if there are more specialists at low latitudes
sp_diet %>% group_by(lat_bin) %>% summarise(mean_diet = mean(gini.index))
# boxplot of specialization
ggplot(sp_diet)+
  geom_boxplot(aes(x=lat_bin, y=gini.index))
# ok they look pretty similar lol
# anova
res.aov <- aov(gini.index ~ lat_bin, data = sp_diet)
summary(res.aov) # it is significant apperently


# Try with habitat data
sp_habitat <- sp_habitat %>% mutate(lat_bin = cut(lat, breaks=c(-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)))
# see if there are more specialists at low latitudes
# boxplot of habitat specialization
ggplot(sp_habitat)+
  geom_boxplot(aes(x=lat_bin, y=Habitat_breadth_IUCN))
# looks like there could be a pattern here
# anova
res.aov <- aov(Habitat_breadth_IUCN ~ lat_bin, data = sp_habitat)
summary(res.aov)
# there does seem to be a significant difference



