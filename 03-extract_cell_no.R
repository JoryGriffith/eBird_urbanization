## This is a script to extract the cell values from the data and add them to the raw data for later

library(terra)
library(sf)
library(tidyverse)
library(auk)


# load global raster script
#GHSL<-rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif") # load raster that is filtered
GHSL<-rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif") # load raster that is filtered
plot(GHSL)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

##### NEED TO ADD AUK UNIQUE AT THIS STEP?? Or maybe just filter out duplicate checklists

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

for (j in 1:length(years)){

  for (i in 1:length(names)){
#  tryCatch(
  dat <- read.delim(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", names[i], "_", years[j], "_filt.txt", sep=""),
               na.strings="", header=TRUE)
  # turn into spatvector
  if(nrow(dat)==0) next
 
  dat <- auk_unique( # remove duplicated group checklists
    dat,
    group_id = "GROUP.IDENTIFIER",
    checklist_id = "SAMPLING.EVENT.IDENTIFIER",
    species_id = "SCIENTIFIC.NAME",
    observer_id = "OBSERVER.ID",
    checklists_only = FALSE
  )
  
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
    rm(dat)
    rm(vect)
    rm(vect2)
  #})
}
  print(paste("finished", years[j]))
}
library(beepr)
beep()
beep_on_error()



