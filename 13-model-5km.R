### Modelling with the scaled up data ######

library(tidyverse)
library(terra)
library(sf)

###########################
# Preparing regular data for modelling
###########################

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

write_csv(datFINAL, "modeling_data_5km.csv")


###################################
## Preparing seasonal data for modelling 
###################################

# load data
summer <- read.csv("summer_richness_summary_5km.csv")
winter <- read.csv("winter_richness_summary_5km.csv")

# extract urbanization values for each
GHSL_5km <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_5km_3levels.tif")
plot(GHSL_5km)

summer$x <- xFromCell(GHSL_5km, summer$cell_5km) # extract the coordinates from the cells
summer$y <- yFromCell(GHSL_5km, summer$cell_5km)

summer$urban <- as.data.frame(terra::extract(GHSL_5km, summer[,c(9:10)]))$SMOD_global
summer$season <- "summer" # add season

# winter
winter$x <- xFromCell(GHSL_5km, winter$cell_5km) # extract the coordinates from the cells
winter$y <- yFromCell(GHSL_5km, winter$cell_5km)

winter$urban <- as.data.frame(terra::extract(GHSL_5km, winter[,c(9:10)]))$SMOD_global
winter$season <- "winter"

season_dat <- rbind(summer, winter) # put them together
# threshold to 102
season_dat_filt <- season_dat %>% filter(number_checklists >=102)
season_dat_filt <- season_dat_filt %>% filter(!is.na(urban))

season_dat_filt %>% group_by(season) %>% summarise(n=n()) # 12,359 summer and 15,254 winter
season_dat_filt %>% group_by(urban) %>% summarise(n=n())

#####################
# extract continent 
continents <- st_read("/Volumes/Backup/eBird/continent-poly/Continents.shp")
#plot(continents)

dat_sf <- st_as_sf(season_dat_filt, coords=c('x', "y"), crs=st_crs(continents))

dat_cont_seas <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature) # joining by nearest feature

###################
# extract biome 
biomes <- st_read("/Volumes/Backup/eBird/wwf_biomes/wwf_terr_ecos.shp")

dat_withbiome_seas <- st_join(dat_cont_seas, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)

# create seperate columns for lat long again
datFINAL_seas <- as.data.frame(dat_withbiome_seas[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
                                                        lat = sf::st_coordinates(.)[,2]))

datFINAL_seas <- datFINAL_seas %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))

# save final data as csv
write_csv(datFINAL_seas, "season_model_data_5km.csv")


#################################
#### Modeling regular data ########
############################

full.dat <- read.csv("modeling_data_5km.csv")
# need to rerun this because it did not save correctly

full.dat$BIOME <- as.factor(full.dat$BIOME)
full.dat$urban <- as.factor(full.dat$urban)

full.dat %>% group_by(urban) %>% summarise(n=n())
summary(full.dat)
hist(full.dat$total_SR, breaks=50)

hist(log(full.dat$total_SR))
hist(sqrt(full.dat$total_SR), breaks=50) # this looks pretty good
hist(full.dat$number_checklists)

# Run models
mod1 <- lm(total_SR ~ abs(lat) * urban + hemisphere + CONTINENT +
             abs(lat):CONTINENT + BIOME + log(number_checklists), dat.full)


dat$abslat <- abs(dat$lat)

mod2 <- lm(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                   BIOME + log(number_checklists), dat.full) # square root transform total species richness

# Plot results of model
predicted <- ggpredict(mod2, terms = c("abslat", "urban")) 

results.plot <-
  plot(predicted, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, line.size=1.5, 
       show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.spacing.y = unit(1, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

## Try to run model with spatial autocorrelation
dat.full.sf <- st_as_sf(dat.full, coords=c('long', "lat"))

dat.full$residuals <- residuals(mod2)
dat.full$fitted <- fitted(mod2)


dat.full.nb <- dnearneigh(dat.full.sf, d1=0, d2=200) # calculate distances
dat.full.lw <- nb2listw(dat.full.nb, style = "W", zero.policy = TRUE) # turn from matrix to list
# Moran's test
lm.morantest(mod2, dat.full.lw, zero.policy = T) # very spatially autocorrelated
beep()
lm.LMtests(mod2, dat.full.lw, test="all", zero.policy = T) # test for spatial error - very significant
beep()

# Try spatially lagged X model
dat.full.slx <- lmSLX(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + 
                        BIOME + log(number_checklists), data = dat.full, listw = dat.full.lw, zero.policy = TRUE)



