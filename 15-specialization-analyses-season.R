############### Specialization analyses for seasonal data ##################
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)

##### Winter #################

winter.model.data <- read.csv("season_modeling_data.csv") %>% filter(season=="winter")

unique(winter.model.data$square)

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c2", "r4c4") # take out r4c1 and r4c4 because they are not in the final modelling data

GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

datalist.years <- list()
datalist.names <- list()

for (i in 1:length(names)){ # come back to r1c4 (4) - no data in r1c4
  for (j in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% winter.model.data$cell) # filter for cells that are in the final dataset
    dat.filt$SCIENTIFIC.NAME<- as.character(dat.filt$SCIENTIFIC.NAME)
    dat.filt$OBSERVATION.DATE<- as.character(dat.filt$OBSERVATION.DATE)
    dat.filt$OBSERVER.ID<- as.character(dat.filt$OBSERVER.ID)
    dat.filt$SAMPLING.EVENT.IDENTIFIER <- as.character(dat.filt$SAMPLING.EVENT.IDENTIFIER)
    dat.filt$OBSERVATION.COUNT <- as.character(dat.filt$OBSERVATION.COUNT)
    dat.filt$GROUP.IDENTIFIER <- as.character(dat.filt$GROUP.IDENTIFIER)
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  wint_dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  wint_dat_uniquesp$x <- xFromCell(GHSL, wint_dat_uniquesp$cell) # extract the coordinates from the cells
  wint_dat_uniquesp$y <- yFromCell(GHSL, wint_dat_uniquesp$cell)
  wint_dat_uniquesp$square <- names[i]
  datalist.names[[i]] <- wint_dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(wint_dat_uniquesp)
}

winter_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(winter_unique_sp$SCIENTIFIC.NAME)) # 8,724 species
write.table(winter_unique_sp, "winter_unique_species.txt", row.names=FALSE)


# Merge with trait data
# Add lat long coordinates to make it easier to bin by latitude
winter_uniquesp <- read.table("winter_unique_species.txt", header=TRUE)
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif")

length(unique(winter_uniquesp$SCIENTIFIC.NAME))


# add in lat long points to more easily bin by latitude
dat_latlong <- st_as_sf(winter_uniquesp, coords=c("x", "y"), crs=st_crs(GHSL))
dat_latlong <- st_transform(dat_latlong, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- as.data.frame(dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                                   lat = sf::st_coordinates(.)[,2]))

# bind this with data in other crs
winter_uniquesp <- cbind(winter_uniquesp, latlong_df[,5:6])

#### Extract urban scores
winter_uniquesp$urban <- as.data.frame(terra::extract(GHSL, winter_uniquesp[,c(3:4)]))$SMOD_global
test <- winter_uniquesp %>% na.omit(urban) # there are no NAs in urban, this is good
# turn urban into 3 categories
winter_uniquesp <- winter_uniquesp %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), "natural", ifelse(urban==30, "urban", "suburban")))

########## Merge with trait data
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")
winter_sp_habitat <- merge(winter_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
length(unique(winter_sp_habitat$SCIENTIFIC.NAME)) # 8,498 species

winter_sp_habitat$abslat <- abs(winter_sp_habitat$lat)
# Try with habitat data
winter_sp_habitat <- winter_sp_habitat %>% mutate(lat_bin = cut(abslat, breaks=abs(c(0, 10, 20, 30, 40, 50, 60, 70, 80))))                                             
# save data with habitat breadth
write.table(winter_sp_habitat, "winter_habitatbreadth.txt", row.names=F)

######## Diet data 
diet<- read.csv("/Volumes/Expansion/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
winter_sp_diet <- merge(winter_uniquesp, diet[, c(9,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific") # merge with species data
length(unique(winter_sp_diet$SCIENTIFIC.NAME)) # 6,904 species

# Bin latitude by 10 degrees
winter_sp_diet$abslat <- abs(winter_sp_diet$lat)
winter_sp_diet <- sp_diet %>% mutate(lat_bin = cut(abslat, breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80)))
# save data with diet
write.table(winter_sp_diet, "winter_dietspec.txt", row.names=F)











############# Summer ##################
summer.model.data <- read.csv("season_modeling_data.csv") %>% filter(season=="summer")


datalist.years <- list()
datalist.names <- list()

for (i in 1:length(names)){ 
  for (j in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% summer.model.data$cell) # filter for cells that are in the final dataset
    dat.filt$SCIENTIFIC.NAME<- as.character(dat.filt$SCIENTIFIC.NAME)
    dat.filt$OBSERVATION.DATE<- as.character(dat.filt$OBSERVATION.DATE)
    dat.filt$OBSERVER.ID<- as.character(dat.filt$OBSERVER.ID)
    dat.filt$SAMPLING.EVENT.IDENTIFIER <- as.character(dat.filt$SAMPLING.EVENT.IDENTIFIER)
    dat.filt$OBSERVATION.COUNT <- as.character(dat.filt$OBSERVATION.COUNT)
    dat.filt$GROUP.IDENTIFIER <- as.character(dat.filt$GROUP.IDENTIFIER)
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  sum_dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  sum_dat_uniquesp$x <- xFromCell(GHSL, sum_dat_uniquesp$cell) # extract the coordinates from the cells
  sum_dat_uniquesp$y <- yFromCell(GHSL, sum_dat_uniquesp$cell)
  sum_dat_uniquesp$square=names[i]
  datalist.names[[i]] <- sum_dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(sum_dat_uniquesp)
}

summer_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(summer_unique_sp$SCIENTIFIC.NAME)) # 10,723 species
write.table(summer_unique_sp, "summer_unique_species.txt", row.names=FALSE)

##########################

# Merge with trait data
# Add lat long coordinates to make it easier to bin by latitude
summer_uniquesp <- read.table("summer_unique_species.txt", header=TRUE)
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif")

length(unique(summer_uniquesp$SCIENTIFIC.NAME))


# add in lat long points to more easily bin by latitude
dat_latlong <- st_as_sf(summer_uniquesp, coords=c("x", "y"), crs=st_crs(GHSL))
dat_latlong <- st_transform(dat_latlong, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- as.data.frame(dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                                   lat = sf::st_coordinates(.)[,2]))

# bind this with data in other crs
summer_uniquesp <- cbind(summer_uniquesp, latlong_df[,5:6])

#### Extract urban scores
summer_uniquesp$urban <- as.data.frame(terra::extract(GHSL, summer_uniquesp[,c(3:4)]))$SMOD_global
test <- summer_uniquesp %>% na.omit(urban) # there are no NAs in urban, this is good
# turn urban into 3 categories
summer_uniquesp <- summer_uniquesp %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), "natural", ifelse(urban==30, "urban", "suburban")))

########## Merge with trait data
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")
summer_sp_habitat <- merge(summer_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
length(unique(summer_sp_habitat$SCIENTIFIC.NAME)) # 8,498 species

summer_sp_habitat$abslat <- abs(summer_sp_habitat$lat)
# Try with habitat data
summer_sp_habitat <- summer_sp_habitat %>% mutate(lat_bin = cut(abslat, breaks=abs(c(0, 10, 20, 30, 40, 50, 60, 70, 80))))                                             
# save data with habitat breadth
write.table(summer_sp_habitat, "summer_habitatbreadth.txt", row.names=F)

######## Diet data 
diet<- read.csv("/Volumes/Expansion/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
summer_sp_diet <- merge(summer_uniquesp, diet[, c(9,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific") # merge with species data
length(unique(summer_sp_diet$SCIENTIFIC.NAME)) # 6,904 species

# Bin latitude by 10 degrees
summer_sp_diet$abslat <- abs(summer_sp_diet$lat)
summer_sp_diet <- sp_diet %>% mutate(lat_bin = cut(abslat, breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80)))
# save data with diet
write.table(summer_sp_diet, "summer_dietspec.txt", row.names=F)







