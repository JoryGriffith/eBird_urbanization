## This is a script to automate the process of summarising the ebird data for each bounding box group
# NOT USING ANYMORE BECAUSE IT SUMMARISES BY YEAR AND WE WANT IT TO SUMMARISE ACROSS YEARS

library(terra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)

# load global raster script
GHSLreproj<-rast("SMOD_global/reprojected.SMOD_global.tif")

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

xmin <- c(-180, -90, 0, 90, 
          -180, -90, -80, -75, -45, 0, 90, 
          -180, -90, 0, 90,
          -180, -90, 0, 90) # xmin for the bounding box

xmax <- c(-90, 0, 90, 180,
          -90, -80, -75, -45, 0, 90, 180,
          -90, 0, 90, 180,
          -90, 0, 90, 180) # xmax for the bounding box

ymin <- c(45, 45, 45, 45,
          0, 0, 0, 0, 0, 0, 0,
          -45, -45, -45, -45,
          -90, -90, -90, -90) # ymin for the bounding box

ymax <- c(90, 90, 90, 90,
          45, 45, 45, 45, 45, 45, 45,
          0, 0, 0, 0, 
          -45, -45, -45, -45) # ymax for the bounding box
bbox <- as.data.frame(cbind(names, xmin, xmax, ymin, ymax))
bbox$xmin <- as.numeric(bbox$xmin)
bbox$xmax <- as.numeric(bbox$xmax)
bbox$ymin <- as.numeric(bbox$ymin)
bbox$ymax <- as.numeric(bbox$ymax)
#write.csv(bbox, "bounding_box_coordinates.csv")


####### 2017
for (i in 19:19){
  
  dat <- read_ebd(paste("eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""))
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  # crop raster by extent
  GHSL_crop <- crop(GHSLreproj, ext(c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSL_crop, xy[,3:4])
  write.table(dat, paste("eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""), row.names=FALSE)
  
  # use lubridate to separate out month
  dat$month <- month(dat$observation_date)
  
  # make dataframe of months and number of checklists
  months <- dat %>% group_by(cell, month) %>% 
    summarise(number_checklists=length(unique(scientific_name)))
  
  months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
  months_wide[is.na(months_wide)]<-0
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat %>%
    group_by(cell) %>%
    summarize(number_checklists=length(unique(sampling_event_identifier)),
              total_SR=length(unique(scientific_name)),
              total_duration=sum(duration_minutes),
              avg_duration=mean(duration_minutes),
              no_months=length(unique(month))
    )
  
  dat_summary$square <- names[i]
  # now need data frame of cell number and urbanization value
  GHSL_values <- as.data.frame(GHSL_crop, xy=TRUE, cells=TRUE)
  
  data <- merge(GHSL_values, dat_summary, by="cell")
  dataframe <- merge(data, months_wide, by="cell")
  
  write.csv(dataframe, paste("eBird_2017_data/custom_bbox/summary/", names[i], "_summary.csv", sep=""))
  
  rm(dat) # remove so it doesnt take up so much memory
}


###### 2018
for (i in 17:19){ # skip R4C1 because there is no data

dat <- read_ebd(paste("eBird_2018_data/custom_bbox/", names[i], "_2018_filt.txt", sep=""))
# turn into spatvector

vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))

# crop raster by extent
GHSL_crop <- crop(GHSLreproj, ext(c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])))

xy=geom(vect)

# get cell number that each point is in
dat$cell=cellFromXY(GHSL_crop, xy[,3:4])
write.table(dat, paste("eBird_2018_data/custom_bbox/", names[i], "_2018_filt.txt", sep=""), row.names=FALSE)

# use lubridate to separate out month
dat$month <- month(dat$observation_date)

# make dataframe of months and number of checklists
months <- dat %>% group_by(cell, month) %>% 
  summarise(number_checklists=length(unique(scientific_name)))

months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
months_wide[is.na(months_wide)]<-0

## aggregate to get number of checklists and raw richness per cell
dat_summary <- dat %>%
  group_by(cell) %>%
  summarize(number_checklists=length(unique(sampling_event_identifier)),
            total_SR=length(unique(scientific_name)),
            total_duration=sum(duration_minutes),
            avg_duration=mean(duration_minutes),
            no_months=length(unique(month))
  )

dat_summary$square <- names[i]
# now need data frame of cell number and urbanization value
GHSL_values <- as.data.frame(GHSL_crop, xy=TRUE, cells=TRUE)

data <- merge(GHSL_values, dat_summary, by="cell")
dataframe <- merge(data, months_wide, by="cell")

write.csv(dataframe, paste("eBird_2018_data/custom_bbox/summary/", names[i], "_summary.csv", sep=""))

rm(dat) # remove so it doesnt take up so much memory
}


######## 2019 come back to r4c2 and r4c3
for (i in 18:19){
  
  dat <- read_ebd(paste("eBird_2019_data/custom_bbox/", names[18], "_2019_filt.txt", sep=""))
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  # crop raster by extent
  GHSL_crop <- crop(GHSLreproj, ext(c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSL_crop, xy[,3:4])
  write.table(dat, paste("eBird_2019_data/custom_bbox/", names[i], "_2019_filt.txt", sep=""), row.names=FALSE)
  
  # use lubridate to separate out month
  dat$month <- month(dat$observation_date)
  
  # make dataframe of months and number of checklists
  months <- dat %>% group_by(cell, month) %>% 
    summarise(number_checklists=length(unique(scientific_name)))
  
  months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
  months_wide[is.na(months_wide)]<-0
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat %>%
    group_by(cell) %>%
    summarize(number_checklists=length(unique(sampling_event_identifier)),
              total_SR=length(unique(scientific_name)),
              total_duration=sum(duration_minutes),
              avg_duration=mean(duration_minutes),
              no_months=length(unique(month))
    )
  
  dat_summary$square <- names[i]
  # now need data frame of cell number and urbanization value
  GHSL_values <- as.data.frame(GHSL_crop, xy=TRUE, cells=TRUE)
  
  data <- merge(GHSL_values, dat_summary, by="cell")
  dataframe <- merge(data, months_wide, by="cell")
  
  write.csv(dataframe, paste("eBird_2019_data/custom_bbox/summary/", names[i], "_summary.csv", sep=""))
  
  rm(dat) # remove so it doesnt take up so much memory
}




####### 2020
for (i in 17:19){
  
  dat <- read_ebd(paste("eBird_2020_data/custom_bbox/", names[i], "_2020_filt.txt", sep=""))
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  # crop raster by extent
  GHSL_crop <- crop(GHSLreproj, ext(c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell=cellFromXY(GHSL_crop, xy[,3:4])
  write.table(dat, paste("eBird_2020_data/custom_bbox/", names[i], "_2020_filt.txt", sep=""), row.names=FALSE)
  
  # use lubridate to separate out month
  dat$month <- month(dat$observation_date)
  
  # make dataframe of months and number of checklists
  months <- dat %>% group_by(cell, month) %>% 
    summarise(number_checklists=length(unique(scientific_name)))
  
  months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
  months_wide[is.na(months_wide)]<-0
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat %>%
    group_by(cell) %>%
    summarize(number_checklists=length(unique(sampling_event_identifier)),
              total_SR=length(unique(scientific_name)),
              total_duration=sum(duration_minutes),
              avg_duration=mean(duration_minutes),
              no_months=length(unique(month))
    )
  
  dat_summary$square <- names[i]
  # now need data frame of cell number and urbanization value
  GHSL_values <- as.data.frame(GHSL_crop, xy=TRUE, cells=TRUE)
  
  data <- merge(GHSL_values, dat_summary, by="cell")
  dataframe <- merge(data, months_wide, by="cell")
  
  write.csv(dataframe, paste("eBird_2020_data/custom_bbox/summary/", names[i], "_summary.csv", sep=""))
  
  rm(dat) # remove so it doesnt take up so much memory
}


####### 2021 # skipped 2, 3, 4
for (i in 3:4){
  
  dat <- read_ebd(paste("eBird_2021_data/custom_bbox/", names[i], "_2021_filt.txt", sep=""))
  # turn into spatvector
  
  #vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
  
  # crop raster by extent
  #GHSL_crop <- crop(GHSLreproj, ext(c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])))
  
  #xy=geom(vect)
  
  # get cell number that each point is in
  #dat$cell=cellFromXY(GHSL_crop, xy[,3:4])
  #write.table(dat, paste("eBird_2021_data/custom_bbox/", names[i], "_2021_filt.txt", sep=""), row.names=FALSE)
  
  # use lubridate to separate out month
  dat$month <- month(dat$observation_date)
  
  # make dataframe of months and number of checklists
  months <- dat %>% group_by(cell, month) %>% 
    summarise(number_checklists=length(unique(scientific_name)))
  
  months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
  months_wide[is.na(months_wide)]<-0
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat %>%
    group_by(cell) %>%
    summarize(number_checklists=length(unique(sampling_event_identifier)),
              total_SR=length(unique(scientific_name)),
              total_duration=sum(duration_minutes),
              avg_duration=mean(duration_minutes),
              no_months=length(unique(month))
    )
  
  dat_summary$square <- names[i]
  # now need data frame of cell number and urbanization value
  GHSL_values <- as.data.frame(GHSL_crop, xy=TRUE, cells=TRUE)
  
  data <- merge(GHSL_values, dat_summary, by="cell")
  dataframe <- merge(data, months_wide, by="cell")
  
  write.csv(dataframe, paste("eBird_2021_data/custom_bbox/summary/", names[i], "_summary.csv", sep=""))
  
  rm(dat) # remove so it doesnt take up so much memory
}


