## This is a script to extract the cell values from the data and add them to the raw data for later

library(terra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)
library(beepr)

# load global raster script
#GHSL<-rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif") # load raster that is filtered
GHSL<-rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif") # load raster that is filtered
plot(GHSL)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

##xmin <- c(-180, -90, 0, 90, 
#          -180, -90, -80, -75, -45, 0, 90, 
#          -180, -90, 0, 90,
#          -180, -90, 0, 90) # xmin for the bounding box
#
#xmax <- c(-90, 0, 90, 180,
#          -90, -80, -75, -45, 0, 90, 180,
#          -90, 0, 90, 180,
#          -90, 0, 90, 180) # xmax for the bounding box
#
#ymin <- c(45, 45, 45, 45,
#          0, 0, 0, 0, 0, 0, 0,
#          -45, -45, -45, -45,
#          -90, -90, -90, -90) # ymin for the bounding box
#
#ymax <- c(90, 90, 90, 90,
#          45, 45, 45, 45, 45, 45, 45,
#          0, 0, 0, 0, 
#          -45, -45, -45, -45) # ymax for the bounding box
#bbox <- as.data.frame(cbind(names, xmin, xmax, ymin, ymax))
#bbox$xmin <- as.numeric(bbox$xmin)
#bbox$xmax <- as.numeric(bbox$xmax)
#bbox$ymin <- as.numeric(bbox$ymin)
#bbox$ymax <- as.numeric(bbox$ymax)
#write.csv(bbox, "bounding_box_coordinates.csv")
#bbox <- read.csv("bounding_box_coordinates.csv")
#dat <- read.delim("/Volumes/Expansion/eBird/eBird_2017_data/custom_bbox/r4c2_2017_filt.txt", header=TRUE, na.strings="")


years <- c(2017, 2018, 2019, 2020, 2021, 2022)

for (j in 1:length(years)){

  for (i in 1:length(names)){
#  tryCatch(
  dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_",years[1], "_data/custom_bbox/", names[1], "_", years[1], "_filt.txt", sep=""), header=TRUE, na.strings="")
  # turn into spatvector
  
  vect <- st_as_sf(dat, crs=st_crs(4326), coords=c("LONGITUDE","LATITUDE"))
  vect2 <- st_transform(vect, crs=crs(GHSL))
  
 
  xy=st_coordinates(vect2)
  
  # get cell number that each point is in
  dat$cell<-cellFromXY(GHSL, xy)
  # also get coordinates for the midpoint of each cell
  write.table(dat, paste("/Volumes/Backup/eBird/eBird_",years[j],"_data/custom_bbox/", names[i], "_",years[j], "_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
  
#  error = function(e){
 #   message(paste("An error occurred for item", i, ":\n"), e)
    
  #})
}
  print(paste("finished", years[j]))
}








