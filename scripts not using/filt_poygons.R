# Filter out all of the data for each of the polygons
library(auk)
library(raster)
library(sf)
library(tidyverse)
library(beepr)

urban_squares <- st_read("GHSL_data_54009_shapefile/GHSL2_0_MWD_L1_tile_schema_land.shp")
# filter out C9 data because already did it
urban_squares2 <- urban_squares <- urban_squares_trans %>% filter(!grepl("_C9",tile_id))
# load ebird 2021 data
eBird_2021 <- auk_ebd("eBird_2021_data/ebd_2021.txt")

names <- paste("eBird_2021_data/polyfilt", urban_squares2$tile_id, "_data_unfilt.txt", sep="") # for saving data filtered by bbox
names2 <- paste("eBird_2021_data/polyfilt", urban_squares2$tile_id, "_data_filt.txt", sep="") # for saving data filtered by columns
names3 <- paste("eBird_2021_data/polyfilt", urban_squares2$tile_id, "_data_polyfilt.csv", sep="") # for saving data filtered by polygon
names4 <- urban_squares$tile_id


# try to load north american data
ebd <- read_ebd("eBird_2021_data/ebd_2021_NAm_filt.txt")
