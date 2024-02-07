### Modelling with the scaled up data ######
library(tidyverse)
library(terra)
library(sf)
library(elevatr)
library(marginaleffects)

###########################
# Preparing regular data for modelling
###########################

## Load summary data
summary <- read.csv("global_richness_summary_5km.csv")

### Use coverage to threshold the data
coverage <- read.csv("coverage_top500_5km.csv")

quantile(coverage$sampsize_95, 0.95) # 103

# Remove data with a sample size less than 102
summary_filt95 <- summary %>% filter(number_checklists >= 103) # went from 569861 to 57873


#### Now I need to load the urbanization data and figure out a threshold for that
GHSL_5km <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_5km_cellsize.tif")

plot(GHSL_5km)



##### Now extract the urbanization score for the summarised values
summary_filt95$x <- xFromCell(GHSL_5km, summary_filt95$cell_5km_scale) # extract the coordinates from the cells
summary_filt95$y <- yFromCell(GHSL_5km, summary_filt95$cell_5km_scale)

summary_filt95$urban <- as.data.frame(terra::extract(GHSL_5km, summary_filt95[,c(20:21)]))$SMOD_global

summary_filt95 %>% group_by(urban) %>% summarise(n=n())

summary_filt <- summary_filt95 %>% filter(!is.na(urban)) # remove NAs - 40366
summary_filt$urban2 <- factor(summary_filt$urban, levels = c("1", "2", "3"),
                               labels = c("Natural", "Suburban", "Urban")) # rename as natural, suburban, and urban

# Turn into sf (spatial) object
summary_filt_sf <- st_as_sf(summary_filt, coords=c("x", "y"), crs=st_crs(GHSL_5km))

## Add elevation
# need to convert the points back to lat long
summary_latlong <- st_transform(summary_filt_sf, crs=st_crs(4326))
latlong_df <- summary_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                         lat = sf::st_coordinates(.)[,2])

latlong_df <- get_elev_point(latlong_df[,c(23,24, 1:21)], prj=crs(summary_latlong), src="aws", overwrite=TRUE)


dat_summary95 <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL_5km)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                                  y = sf::st_coordinates(.)[,2]))


######## Extract continent
#continents <- st_read("/Volumes/Backup/eBird/continent-poly/Continents.shp")
##plot(continents)
#
#dat_sf <- st_as_sf(dat_summary95, coords=c('long', "lat"), crs=st_crs(continents))
#
#dat_cont <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature)
#
#
######### Extract biome
#biomes <- st_read("/Volumes/Expansion/eBird/wwf_biomes/wwf_terr_ecos.shp")
#class(biomes) # sf and data frame
#
#dat_withbiome <- st_join(dat_sf, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)
#
#
## create seperate columns for lat long again
#datFINAL <- as.data.frame(dat_withbiome[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
 #                                                       lat = sf::st_coordinates(.)[,2]))

#summary(datFINAL)
# save as csv
#write_csv(datFINAL, "modeling_data_5km.csv")

################
### Add precipitation
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
dat_summary95$precip <- as.data.frame(terra::extract(precip, dat_summary95[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12



######## Assign Hemisphere
dat_summary95$hemisphere <- "northern"
dat_summary95$hemisphere[dat_summary95$lat<0]<-"southern"


dat_summary95$abslat <- abs(dat_summary95$lat)

write_csv(dat_summary95, "modeling_data_5km.csv")





########################### MODELLING ###################
full.dat <- read.csv("modeling_data_5km.csv")


full.dat %>% group_by(urban) %>% summarise(n=n())
summary(full.dat)
hist(full.dat$total_SR, breaks=50)

hist(log(full.dat$total_SR))
hist(sqrt(full.dat$total_SR), breaks=50) # this looks pretty good
hist(full.dat$number_checklists)

# filter out latitudes without urban cells
full.dat <- full.dat %>% filter(lat <= 70 & lat >=-55)


ggplot(full.dat, aes(x=abslat, y=total_SR, color=urban2))+
#  geom_point(alpha=0.1)+
  geom_smooth(method="lm")


# something is going wrong with the elevation
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
             precip + log(number_checklists) + elevation, full.dat)

# Plot results of model
square <- function(x){
  x^2
} 

## plot results
predicted.5km<-avg_predictions(mod1, by=c("abslat", "urban2"), transform=square, 
                              newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                 urban2=c("Natural", "Suburban", "Urban")))


mod_5km.plot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(full.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.5km, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.5km, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mod_5km.plot

plot_slopes(mod1, variables="abslat", condition=c("urban2"))

























###################################
## Preparing seasonal data for modelling 
###################################

# load data
summer <- read.csv("summer_richness_summary_5km.csv")
winter <- read.csv("winter_richness_summary_5km.csv")
# ok something went wrong here


# extract urbanization values for each
GHSL_5km <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_5km_cellsize.tif")
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
season_dat_filt <- season_dat %>% filter(number_checklists >=103)
season_dat_filt <- season_dat_filt %>% filter(!is.na(urban))

season_dat_filt %>% group_by(season) %>% summarise(n=n()) # 12,359 summer and 15,254 winter
season_dat_filt %>% group_by(urban) %>% summarise(n=n())


season_dat_filt$urban2 <- factor(season_dat_filt$urban, levels = c("1", "2", "3"),
                              labels = c("Natural", "Suburban", "Urban")) # rename as natural, suburban, and urban

#####################
# extract continent 
#continents <- st_read("/Volumes/Backup/eBird/continent-poly/Continents.shp")
##plot(continents)
#
#dat_sf <- st_as_sf(season_dat_filt, coords=c('x', "y"), crs=st_crs(continents))
#
#dat_cont_seas <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature) # joining by nearest feature

###################
# extract biome 
#biomes <- st_read("/Volumes/Backup/eBird/wwf_biomes/wwf_terr_ecos.shp")
#
#dat_withbiome_seas <- st_join(dat_cont_seas, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)
#
## create seperate columns for lat long again
#datFINAL_seas <- as.data.frame(dat_withbiome_seas[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
#                                                        lat = sf::st_coordinates(.)[,2]))
#
#datFINAL_seas <- datFINAL_seas %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))
#
## save final data as csv
#write_csv(datFINAL_seas, "season_model_data_5km.csv")

################
### Add elevation
# Turn into sf (spatial) object
season_dat_filt_sf <- st_as_sf(season_dat_filt, coords=c("x", "y"), crs=st_crs(GHSL_5km))

summary_latlong <- st_transform(season_dat_filt_sf, crs=st_crs(4326)) # change crs
latlong_df <- summary_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                         lat = sf::st_coordinates(.)[,2])

latlong_df <- get_elev_point(latlong_df[,c(13,14,2:11)], prj=crs(summary_latlong), src="aws", overwrite=TRUE)


dat_season <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL_5km)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                                      y = sf::st_coordinates(.)[,2]))



### Add precipitation
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
dat_season$precip <- as.data.frame(terra::extract(precip, dat_season[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12



######## Assign Hemisphere
dat_season$hemisphere <- "northern"
dat_season$hemisphere[dat_season$lat<0]<-"southern"


dat_season$abslat <- abs(dat_season$lat)


write_csv(dat_season, "season_model_data_5km.csv")


#################################
#### Modeling seasonal data ########
############################

season.dat <- read.csv("season_model_data_5km.csv")
# need to rerun this because it did not save correctly

season.dat %>% group_by(urban) %>% summarise(n=n())
summary(season.dat)
hist(season.dat$total_SR, breaks=50)

hist(log(season.dat$total_SR))
hist(sqrt(season.dat$total_SR), breaks=50) # this looks pretty good
hist(season.dat$number_checklists)


# Run models
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2 * season + hemisphere + 
             precip + log(number_checklists) + elevation, season.dat)
# added hemisphere as a intercept effect because model was rank deficient


square <- function(x){
  x^2
} 

## plot results
predicted.season.5km<-avg_predictions(mod1, by=c("abslat", "urban2", "season"), transform=square, 
                                      newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                         season = c("summer", "winter")))


mod_5km.plot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(season.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season.5km, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season.5km, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())+
  facet_wrap(~season)
mod_5km.plot

plot_slopes(mod1, variables="abslat", condition=c("urban2", "season"))







######### Old code for categorizing by urbanization


# aggregate
# make 2 layers, one with the percentage of water (NA) cells and one with the mean values
# calculate percentage of NAs
#GHSL_percNA <- aggregate(is.na(GHSLreproj$SMOD_global), fact=5, mean) * 100
#plot(GHSL_percNA)
## this looks pretty good
#
## now do another layer that is the mean of the rest
#GHSL_5km <- aggregate(GHSLreproj, fact=5, mean, na.rm=TRUE)
##GHSL_5km.test <- GHSL_5km
## remove the cells that are more than 20% water
#GHSL_5km[(GHSL_percNA > 20)] <- NA
#
#summary(values(GHSL_percNA))
#summary(values(GHSL_5km))
##summary(values(GHSL_5km.test)) # ok that worked, removed squares that have over 20% water
#
## Now I need to categorize the rest of the data
## I think I should also do it by percentages
## percentage of category 1 (rural)
#GHSL_perc1 <- aggregate(GHSLreproj$SMOD_global==1, fact=5, mean) * 100
#plot(GHSL_perc1)
#
#GHSL_perc2 <- aggregate(GHSLreproj$SMOD_global==2, fact=5, mean) * 100
#plot(GHSL_perc2)
#
#GHSL_perc3 <- aggregate(GHSLreproj$SMOD_global==3, fact=5, mean) * 100
#plot(GHSL_perc3)
#
## Categorize cells that are more than 60% of one thing as that thing
#GHSL_5km[(GHSL_perc1 > 60)] <- 1 # categorize greater than 60% rural as rural
#GHSL_5km[(GHSL_perc2 > 60)] <- 2 # categorize greater than 60% peri-urban as peri-urban
#GHSL_5km[(GHSL_perc3 > 60)] <- 3 # categorize greater than 60% urban as urban
## categorize everything else as NA
#GHSL_5km[(GHSL_5km$SMOD_global>1) & (GHSL_5km$SMOD_global<2)] <- NA
#GHSL_5km[(GHSL_5km$SMOD_global>2) & (GHSL_5km$SMOD_global<3)] <- NA
#unique(values(GHSL_5km))
#summary(values(GHSL_5km))
#
#plot(GHSL_5km)
##GHSL.test <- GHSL_5km
##GHSL.test[(GHSL.test$SMOD_global<=1.2)] <- 1
##GHSL.test[(GHSL.test$SMOD_global>1.2) & (GHSL.test$SMOD_global<=1.9)] <- NA
##GHSL.test[(GHSL.test$SMOD_global>=1.9) & (GHSL.test$SMOD_global<=2.1)] <- 2
##GHSL.test[(GHSL.test$SMOD_global>2.1) & (GHSL.test$SMOD_global<2.8)] <- NA
##GHSL.test[(GHSL.test$SMOD_global>=2.8)] <- 3
#
#writeRaster(GHSL_5km, "/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_3levels.tif")
#


