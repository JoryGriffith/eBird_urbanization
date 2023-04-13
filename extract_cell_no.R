## This is a script to extract the cell values from the data and add them to the raw data for later

library(terra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)
library(beepr)

# load global raster script
GHSLreproj<-rast("SMOD_global/reprojected.SMOD_global.tif")

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


####### 2017
for (i in 17:17){
  
  dat <- read.table(paste("eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""), header=TRUE)
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSLreproj, xy[,3:4])
  # also get coordinates for the midpoint of each cell
  write.table(dat, paste("eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
}


###### 2018
for (i in 17:19){
  
  dat <- read.table(paste("eBird_2018_data/custom_bbox/", names[i], "_2018_filt.txt", sep=""), header=TRUE)
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSLreproj, xy[,3:4])
  write.table(dat, paste("eBird_2018_data/custom_bbox/", names[i], "_2018_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
}


######## 2019 (skipped r4c3)
for (i in 18:19){
  
  dat <- read.table(paste("eBird_2019_data/custom_bbox/", names[i], "_2019_filt.txt", sep=""), header=TRUE)
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSLreproj, xy[,3:4])
  write.table(dat, paste("eBird_2019_data/custom_bbox/", names[i], "_2019_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
}


####### 2020
for (i in 17:19){
  
  dat <- read.table(paste("eBird_2020_data/custom_bbox/", names[i], "_2020_filt.txt", sep=""), header=TRUE)
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSLreproj, xy[,3:4])
  write.table(dat, paste("eBird_2020_data/custom_bbox/", names[i], "_2020_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
}


####### 2021
for (i in 1:16){
    
    dat <- read.table(paste("eBird_2021_data/custom_bbox/", names[i], "_2021_filt.txt", sep=""), header=TRUE)
    # turn into spatvector
    
    vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
    
    xy=geom(vect)
    
    # get cell number that each point is in
    dat$cell=cellFromXY(GHSLreproj, xy[,3:4])
    write.table(dat, paste("eBird_2021_data/custom_bbox/", names[i], "_2021_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
  }


