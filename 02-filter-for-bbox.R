##### This is to filter for custom bboxes for the different years so I can load the data into R

library(terra)
library(sf)
library(auk)
library(beepr)

# load datasets
eBird_2017 <- auk_ebd("/Volumes/Expansion/eBird/eBird_2017_data/ebd_2017.txt")
eBird_2018 <- auk_ebd("/Volumes/Expansion/eBird/eBird_2018_data/ebd_2018.txt")
eBird_2019 <- auk_ebd("/Volumes/Expansion/eBird/eBird_2019_data/ebd_2019.txt")
eBird_2020 <- auk_ebd("/Volumes/Expansion/eBird/eBird_2020_data/ebd_2020.txt")
eBird_2021 <- auk_ebd("/Volumes/Expansion/eBird/eBird_2021_data/ebd_2021.txt")
eBird_2022 <- auk_ebd("/Volumes/Expansion/eBird/eBird_2022_data/ebd_2022.txt")

# load columns I want
cols <- c("latitude", "longitude", "group identifier", "sampling event identifier",
          "scientific name", "observation count", "observer_id", "observation_date", "duration_minutes")

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4") # names of the bounding boxes

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
# filter data for the custom bbox using a loop (for 2018)

# 2017
for (i in 1:19){
  
  eBird_2017 %>% 
    auk_bbox(bbox=c(bbox$xmin[i], bbox$ymin[i], bbox$xmax[i], bbox$ymax[i])) %>% # filter out square
    auk_filter(file = paste("/Volumes/Expansion/eBird/eBird_2017_data/custom_bbox/", names[i], "_2017_unfilt.txt", sep="")) # make dataframe
  
  # then filter out columns I want
  
  auk_ebd(paste("/Volumes/Expansion/eBird/eBird_2017_data/custom_bbox/", names[i], "_2017_unfilt.txt", sep="")) %>% 
    auk_select(select = cols, file = paste("/Volumes/Expansion/eBird/eBird_2017_data/custom_bbox/", names[i], "_2017_filt.txt", sep=""))
  
}

# 2018
for (i in 1:19){
  
  eBird_2018 %>% 
    auk_bbox(bbox=c(bbox$xmin[i], bbox$ymin[i], bbox$xmax[i], bbox$ymax[i])) %>% # filter out square
    auk_filter(file = paste("/Volumes/Expansion/eBird/eBird_2018_data/custom_bbox/", names[i], "_2018_unfilt.txt", sep="")) # make dataframe
  
  # then filter out columns I want
  
  auk_ebd(paste("/Volumes/Expansion/eBird/eBird_2018_data/custom_bbox/", names[i], "_2018_unfilt.txt", sep="")) %>% 
    auk_select(select = cols, file = paste("/Volumes/Expansion/eBird/eBird_2018_data/custom_bbox/", names[i], "_2018_filt.txt", sep=""))
  
}

# 2019
for (i in 1:19){
  
  eBird_2019 %>% 
    auk_bbox(bbox=c(bbox$xmin[i], bbox$ymin[i], bbox$xmax[i], bbox$ymax[i])) %>% # filter out square
    auk_filter(file = paste("/Volumes/Expansion/eBird/eBird_2019_data/custom_bbox/", names[i], "_2019_unfilt.txt", sep="")) # make dataframe
  
  # then filter out columns I want
  
  auk_ebd(paste("/Volumes/Expansion/eBird/eBird_2019_data/custom_bbox/", names[i], "_2019_unfilt.txt", sep="")) %>% 
    auk_select(select = cols, file = paste("/Volumes/Expansion/eBird/eBird_2019_data/custom_bbox/", names[i], "_2019_filt.txt", sep=""))
  
}


# 2020
for (i in 1:19){
  
  eBird_2020 %>% 
    auk_bbox(bbox=c(bbox$xmin[i], bbox$ymin[i], bbox$xmax[i], bbox$ymax[i])) %>% # filter out square
    auk_filter(file = paste("/Volumes/Expansion/eBird/eBird_2020_data/custom_bbox/", names[i], "_2020_unfilt.txt", sep="")) # make dataframe
  
  # then filter out columns I want
  
  auk_ebd(paste("/Volumes/Expansion/eBird/eBird_2020_data/custom_bbox/", names[i], "_2020_unfilt.txt", sep="")) %>% 
    auk_select(select = cols, file = paste("/Volumes/Expansion/eBird/eBird_2020_data/custom_bbox/", names[i], "_2020_filt.txt", sep=""))
  
}


# 2021
for (i in 1:19){
  
  eBird_2021 %>% 
    auk_bbox(bbox=c(bbox$xmin[i], bbox$ymin[i], bbox$xmax[i], bbox$ymax[i])) %>% # filter out square
    auk_filter(file = paste("/Volumes/Expansion/eBird/eBird_2021_data/custom_bbox/", names[i], "_2021_unfilt.txt", sep="")) # make dataframe
  
  # then filter out columns I want
  
  auk_ebd(paste("/Volumes/Expansion/eBird/eBird_2021_data/custom_bbox/", names[i], "_2021_unfilt.txt", sep="")) %>% 
    auk_select(select = cols, file = paste("/Volumes/Expansion/eBird/eBird_2021_data/custom_bbox/", names[i], "_2021_filt.txt", sep=""))
  
}

# 2022
for (i in 1:19){
  
  eBird_2022 %>% 
    auk_bbox(bbox=c(bbox$xmin[i], bbox$ymin[i], bbox$xmax[i], bbox$ymax[i])) %>% # filter out square
    auk_filter(file = paste("/Volumes/Expansion/eBird/eBird_2022_data/custom_bbox/", names[i], "_2022_unfilt.txt", sep="")) # make dataframe
  
  # then filter out columns I want
  
  auk_ebd(paste("/Volumes/Expansion/eBird/eBird_2022_data/custom_bbox/", names[i], "_2022_unfilt.txt", sep="")) %>% 
    auk_select(select = cols, file = paste("/Volumes/Expansion/eBird/eBird_2022_data/custom_bbox/", names[i], "_2022_filt.txt", sep=""))
  
}
