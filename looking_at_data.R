library(tidyverse)
library(terra)

#### Looking at why there is no data in south africa
i <- 2021
dat <- read.table("/Volumes/Expansion/eBird/eBird_2021_data/custom_bbox/R3C3_2021_filt.txt", header=TRUE)

summary <- dat %>% group_by(cell) %>% summarise(n=n(),
                                                number_checklists=length(unique(sampling_event_identifier)),
                                                lat=mean(latitude),
                                                long=mean(longitude)) 



# load GHSL raster
GHSL_filt <- rast("/Volumes/Expansion/eBird/SMOD_global/urbanization.tif")
# extract cell numbers


# load summary filtered data
summary_filt <- read.csv("5yr_summary/summary_thresholded.csv")

check <- summary_filt %>% filter(cell==417884558)

vect <- vect(dat, crs=crs(GHSL_filt),geom=c("longitude","latitude"))

xy=geom(vect)

# get cell number that each point is in
dat$ext <- terra::extract(GHSL_filt, xy[,3:4])

# try again with summary data
vect <- vect(summary, crs=crs(GHSL_filt),geom=c("long","lat"))

xy=geom(vect)

# get cell number that each point is in
summary$ext <- terra::extract(GHSL_filt, xy[,3:4])

# There are places in south africa that should have points which means that something is going on...need to look more into this
# maybe need to rerun code and check if there is something happening with the cell extraction or summarising, not sure




