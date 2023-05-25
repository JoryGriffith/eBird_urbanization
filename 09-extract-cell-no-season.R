#### Here, I will extract the cell number for the seasonal data
library(terra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)
library(beepr)

# load global raster script
#GHSLreproj<-rast("/Volumes/Expansion/eBird/SMOD_global/reprojected.SMOD_global.tif")
GHSLreproj<-rast("/Volumes/Backup/eBird/SMOD_global/reprojected.SMOD_global.tif")

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

# Summer
for (j in 1:length(years)){
  for (i in 1:length(names)){
   
     dat <- read.delim(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), header=TRUE, na.strings="")
    # turn into spatvector
    vect <- vect(dat, crs=crs(GHSLreproj),geom=c("LONGITUDE","LATITUDE"))
    
    xy=geom(vect)
    
    # get cell number that each point is in
    dat$cell<-cellFromXY(GHSLreproj, xy[,3:4])
    # also get coordinates for the midpoint of each cell
    write.table(dat, paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
    
  }
  print(paste("finished", years[j]))
}


# Winter
for (j in 1:length(years)){
  for (i in 1:length(names)){
    
    dat <- read.delim(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), header=TRUE, na.strings="")
    # turn into spatvector
    vect <- vect(dat, crs=crs(GHSLreproj),geom=c("LONGITUDE","LATITUDE"))
    
    xy=geom(vect)
    
    # get cell number that each point is in
    dat$cell<-cellFromXY(GHSLreproj, xy[,3:4])
    # also get coordinates for the midpoint of each cell
    write.table(dat, paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
    
  }
  print(paste("finished", years[j]))
}



