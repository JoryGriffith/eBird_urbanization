## In this script, I will start running some models of how the latitudinal diversity gradient changes with urbanization
library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world)

# load thresholded summary data
dat <- read.csv("5yr_summary/summary_thresholded.csv")
dat <- rename(dat, lat = y, long = x)

# First I need to get the other data that I want to include in the model


######## Assign hemisphere
dat <- dat %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))

######### Assign continent
dat_sf <- st_as_sf(dat, coords=c('long', "lat"), crs=st_crs(world))

joined <- st_join(dat_sf, world)

dat_joined <- as.data.frame(joined[,c(1:22,41)] %>% mutate(long = sf::st_coordinates(.)[,1],
                                                           lat = sf::st_coordinates(.)[,2]))# just keep country name
# extract continent using country name
dat_joined$continent <- countrycode(sourcevar = dat_joined[,"name_long"],
                                     origin = "country.name",
                                     destination = "continent")


#########################
# Extract biome
# Classifying points into biomes using a terrestrial biomes shapefile
biomes <- st_read("/Volumes/Expansion/eBird/wwf_biomes/wwf_terr_ecos.shp")
class(biomes) # sf and data frame
# look at how many biomes there are
length(unique(biomes$REALM)) # 9 of these
length(unique(biomes$BIOME)) # 16 biomes
length(unique(biomes$ECO_NAME)) # 827 ecoregion names

# want to extract biomes
dat_sf <- st_as_sf(dat_joined, crs=st_crs(biomes))
dat_withbiome <- st_join(dat_sf, biomes[,"BIOME"])

# create seperate columns for lat long again
datFINAL <- dat_withbiome %>% mutate(long = sf::st_coordinates(.)[,1],
                                             lat = sf::st_coordinates(.)[,2])
# have a bunch of NA values for continent and biome, not sure why need to look into that

# save as csv
write.csv(datFINAL, "modeling_data.csv")




