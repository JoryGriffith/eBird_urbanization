########### Script for analyzing whether specialist species are being disproportionately lost at lower latitudes
library(tidyverse)
library(terra)
## First I want to thin the data so that for each cell, there is only one row for each species 

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

model.data <- read.csv("modeling_data.csv") 
# also going to load this because I want to filter for cells that are in the final dataset to save space

GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

datalist.years <- list()
datalist.names <- list()

for (i in 2:length(names)){
  for (j in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", names[i], "_", years[j], "_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% model.data$cell) # filter for cells that are in the final dataset
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  dat_uniquesp$x <- xFromCell(GHSL, dat_uniquesp$cell) # extract the coordinates from the cells
  dat_uniquesp$y <- yFromCell(GHSL, dat_uniquesp$cell)
  datalist.names[[i]] <- dat_uniquesp
}

global_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together

write.csv(global_unique_sp, paste("global_unique_species.csv", sep=""))

### Put together with 






