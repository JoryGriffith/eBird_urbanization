############### Specialization analyses for seasonal data ##################
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)

##### Winter #################

winter.model.data <- read.csv("season_model_data.csv") %>% filter(season=="winter")

unique(winter.model.data$square)

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c2", "r4c4") # take out r4c1 and r4c4 because they are not in the final modelling data

GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

datalist.years <- list()
datalist.names <- list()

for (i in 1:length(names)){ 
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
  wint_dat_uniquesp$square=names[i]
  datalist.names[[i]] <- wint_dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(wint_dat_uniquesp)
}

winter_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(winter_unique_sp$SCIENTIFIC.NAME)) # 10,723 species
write.table(winter_unique_sp, "winter_unique_species.txt", row.names=FALSE)





############# Summer ##################
summer.model.data <- read.csv("season_model_data.csv") %>% filter(season=="summer")


datalist.years <- list()
datalist.names <- list()

for (i in 1:length(names)){ 
  for (j in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), 
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








