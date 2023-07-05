### Modelling with the scaled up data ######

library(tidyverse)
library(terra)

## Load summary data
summary <- read.csv("global_richness_summary_5km.csv")

### Use coverage to threshold the data
coverage <- read.csv("coverage_top500_5km.csv")

quantile(coverage$sampsize_95, 0.95) # 102

# Remove data with a sample size less than 102
summary_filt95 <- summary %>% filter(number_checklists >= 102) # went from 569861 to 57873


#### Now I need to load the urbanization data and figure out a threshold for that
GHSL_5km <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_cellsize.tif")

summary(values(GHSL_5km))
# right now it is an average, so I need to figure out thresholds for each category

plot(GHSL_5km)

# Load original data
GHSLreproj<-rast("/Volumes/Expansion/eBird/SMOD_global/reprojected.SMOD_global.tif")
unique(values(GHSLreproj))
summary(values(GHSLreproj)) # NAs 4029
# recategorize to 1, 2, and 3

GHSLreproj[(GHSLreproj==10)] <- NA
GHSLreproj <- subst(GHSLreproj, c(11,12, 13), 1)
GHSLreproj <- subst(GHSLreproj, c(21,22,23), 2)
GHSLreproj <- subst(GHSLreproj, 30, 3)

# aggregate
# make 2 layers, one with the percentage of water (NA) cells and one with the mean values
# calculate percentage of NAs
GHSL_percNA <- aggregate(is.na(GHSLreproj$SMOD_global), fact=5, mean) * 100
plot(GHSL_percNA)
# this looks pretty good

# now do another layer that is the mean of the rest
GHSL_5km <- aggregate(GHSLreproj, fact=5, mean, na.rm=TRUE)
#GHSL_5km.test <- GHSL_5km
# remove the cells that are more than 20% water
GHSL_5km[(GHSL_percNA > 20)] <- NA

summary(values(GHSL_percNA))
summary(values(GHSL_5km))
#summary(values(GHSL_5km.test)) # ok that worked, removed squares that have over 20% water

# Now I need to categorize the rest of the data
# I think I should also do it by percentages
# percentage of category 1 (rural)
GHSL_perc1 <- aggregate(GHSLreproj$SMOD_global==1, fact=5, mean) * 100
plot(GHSL_perc1)

GHSL_perc2 <- aggregate(GHSLreproj$SMOD_global==2, fact=5, mean) * 100
plot(GHSL_perc2)

GHSL_perc3 <- aggregate(GHSLreproj$SMOD_global==3, fact=5, mean) * 100
plot(GHSL_perc3)

# Categorize cells that are more than 60% of one thing as that thing
GHSL_5km[(GHSL_perc1 > 60)] <- 1 # categorize greater than 60% rural as rural
GHSL_5km[(GHSL_perc2 > 60)] <- 2 # categorize greater than 60% peri-urban as peri-urban
GHSL_5km[(GHSL_perc3 > 60)] <- 3 # categorize greater than 60% urban as urban
# categorize everything else as NA
GHSL_5km[(GHSL_5km$SMOD_global>1) & (GHSL_5km$SMOD_global<2)] <- NA
GHSL_5km[(GHSL_5km$SMOD_global>2) & (GHSL_5km$SMOD_global<3)] <- NA
unique(values(GHSL_5km))
summary(values(GHSL_5km))

plot(GHSL_5km)
#GHSL.test <- GHSL_5km
#GHSL.test[(GHSL.test$SMOD_global<=1.2)] <- 1
#GHSL.test[(GHSL.test$SMOD_global>1.2) & (GHSL.test$SMOD_global<=1.9)] <- NA
#GHSL.test[(GHSL.test$SMOD_global>=1.9) & (GHSL.test$SMOD_global<=2.1)] <- 2
#GHSL.test[(GHSL.test$SMOD_global>2.1) & (GHSL.test$SMOD_global<2.8)] <- NA
#GHSL.test[(GHSL.test$SMOD_global>=2.8)] <- 3

writeRaster(GHSL_5km, "/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_3levels.tif")

##### Now extract the urbanization score for the summarised values
summary_filt95$x <- xFromCell(GHSL_5km, summary_filt95$cell_5km_scale) # extract the coordinates from the cells
summary_filt95$y <- yFromCell(GHSL_5km, summary_filt95$cell_5km_scale)

summary_filt95$urban <- as.data.frame(terra::extract(GHSL_5km, summary_filt95[,c(21:22)]))$SMOD_global

summary_filt95 %>% group_by(urban) %>% summarise(n=n())

summary_filt <- summary_filt95 %>% filter(!is.na(urban)) # remove NAs - 44340

summary_filt <- rename(summary_filt, lat = y, long = x) # rename to lat and long

######## Assign Hemisphere
summary_filt <- summary_filt %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))

######## Extract continent
continents <- st_read("/Volumes/Expansion/eBird/continent-poly/Continents.shp")
#plot(continents)

dat_sf <- st_as_sf(summary_filt, coords=c('long', "lat"), crs=st_crs(continents))

dat_cont <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature)


######## Extract biome
biomes <- st_read("/Volumes/Expansion/eBird/wwf_biomes/wwf_terr_ecos.shp")
class(biomes) # sf and data frame

dat_withbiome <- st_join(dat_cont, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)


# create seperate columns for lat long again
datFINAL <- as.data.frame(dat_withbiome[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
                                                        lat = sf::st_coordinates(.)[,2]))

summary(datFINAL)
# save as csv
write.csv(datFINAL, "modeling_data_5km.csv", row.names=FALSE)






