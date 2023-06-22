############### Code to scale up the analysis to 5x5 km rasters as a sensitivity analysis #####################

library(terra)
library(tidyverse)
library(lubridate)
library(beepr)
library(sf)

# load original reprojected raster file

GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/reprojected.SMOD_global.tif")

GHSL_5km <- aggregate(GHSL, fact=5, fun="mean", cores=4, filename="/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_cellsize.tif")

plot(GHSL_5km)

# yay! created the new raster

# Now I need to aggregate the eBird data into that raster and summarise

########## 1) extract cell no. for each point ########
names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

for (j in 1:length(years)) {

  for (i in 1:length(names)){
  dat <- read.delim(paste("/Volumes/Expansion/eBird/eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""), header=TRUE, na.strings="")
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSL_5km),geom=c("LONGITUDE","LATITUDE"))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell_5km_scale<-cellFromXY(GHSL_5km, xy[,3:4])
  # also get coordinates for the midpoint of each cell
  write.table(dat, paste("/Volumes/Expansion/eBird/eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
}
print(paste("finished", years[j]))
}


############# 2) Summarise by cell no. ########
for (j in 1:length(names)){
  datalist = vector("list", length = length(years))
  # loop for each year
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_", years[i], "_data/custom_bbox/", names[14], "_", years[i], "_filt.txt", sep=""), 
                      header=TRUE) # load data
    
    dat$SCIENTIFIC.NAME <- as.character(dat$SCIENTIFIC.NAME)
    dat$OBSERVATION.DATE <- as.character(dat$OBSERVATION.DATE)
    dat$OBSERVER.ID <- as.character(dat$OBSERVER.ID)
    dat$SAMPLING.EVENT.IDENTIFIER <- as.character(dat$SAMPLING.EVENT.IDENTIFIER)
    dat$OBSERVATION.COUNT <- as.character(dat$OBSERVATION.COUNT)
    dat$GROUP.IDENTIFIER <- as.character(dat$GROUP.IDENTIFIER)
    datalist[[i]] <- dat # put in a list
  }
  # bind lists together
  
  dat2 <- dplyr::bind_rows(datalist)
  
  # summarise
  dat2$month <- month(dat2$OBSERVATION.DATE)
  
  # make dataframe of months and number of checklists
  months <- dat2 %>% group_by(cell_5km, month) %>% 
    summarise(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)))
  
  months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
  months_wide[is.na(months_wide)]<-0
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat2 %>%
    group_by(cell_5km) %>%
    summarize(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)),
              total_SR=length(unique(SCIENTIFIC.NAME)),
              total_duration=sum(DURATION.MINUTES),
              avg_duration=mean(DURATION.MINUTES),
              no_months=length(unique(month))
    )
  
  dat_summary2 <- merge(dat_summary, months_wide, by="cell_5km")
  dat_summary2$square=names[j]
  # save as csv
  write.csv(dat_summary2, paste("5yr_summary_5km/", names[j], "_SR.csv", sep=""))
  print(paste("finished", names[j]))
}

##### put all together
list_csv_files <- list.files(path = "5yr_summary_5km/", pattern="*.csv")

names <- tolower(gsub('_SR.csv', "", list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary_5km/", list_csv_files[i], sep="")))
}

dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)
#save all the summaries as a csv
write.csv(dat, "global_richness_summary_5km.csv", row.names=FALSE)

# find top 500 cells
top_cells <- dat %>% slice_max(total_SR, n=500)
write.csv(top_cells, "top_500_cells_5km.csv", row.names=FALSE)


############# 3) Calculate threshold for inclusion ##############
top_cells_5km <- read.csv("top_500_cells_5km.csv") 



