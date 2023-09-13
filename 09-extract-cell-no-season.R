#### Here, I will extract the cell number for the seasonal data
library(terra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)
library(beepr)

# load global raster script
#GHSLreproj<-rast("/Volumes/Expansion/eBird/SMOD_global/reprojected.SMOD_global.tif")
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# Summer
for (j in 1:length(years)){
  for (i in 1:length(names)){
   
     dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), header=TRUE, na.strings="")
    # turn into spatvector
     
     vect <- st_as_sf(dat, crs=st_crs(4326), coords=c("LONGITUDE","LATITUDE"))
     vect2 <- st_transform(vect, crs=crs(GHSL))
    
     xy=st_coordinates(vect2)
    
    # get cell number that each point is in
    dat$cell<-cellFromXY(GHSL, xy)
    # also get coordinates for the midpoint of each cell
    write.table(dat, paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
    rm(dat)
    rm(vect)
    rm(vect2)
  }
  print(paste("finished", years[j]))
}


# Winter
for (j in 1:length(years)){
  for (i in 1:length(names)){
    
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), header=TRUE, na.strings="")
    # turn into spatvector
    vect <- st_as_sf(dat, crs=st_crs(4326), coords=c("LONGITUDE","LATITUDE"))
    vect2 <- st_transform(vect, crs=crs(GHSL))
    
    xy=st_coordinates(vect2)
    
    # get cell number that each point is in
    dat$cell<-cellFromXY(GHSL, xy)
    # also get coordinates for the midpoint of each cell
    write.table(dat, paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
    rm(dat)
    rm(vect)
    rm(vect2)
  }
  print(paste("finished", years[j]))
}



