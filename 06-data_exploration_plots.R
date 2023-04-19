############# Make plots of the data that has been thresholded and filtered ###########
library(terra)
library(tidyverse)
library(sf)
library(ggpubr)
library(tidyterra)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
#####################################
# Data exploration of the filtered data
# Load csv
#summary_filt <- read.csv("5yr_summary/summary_thresholded.csv")

# try to rasterize to plot
# try to aggregate to make more visible
#plot(thresholded.raster)
#raster.aggregated <- aggregate(thresholded.raster, fact = 10, fun="sum")
#raster.aggregated[is.na(raster.aggregated)] <- 10000
#plot(raster.aggregated)


# plot of urbanization and species richness
ggplot(summary_filt, aes(x=urban, y=total_SR))+
  geom_point()+
  geom_smooth(method="lm")
# definitely going down which is great

# plot of species richness and number of checklists
ggplot(summary_filt, aes(x=log(number_checklists), y=log(total_SR)))+
  geom_point()+
  geom_smooth(method="lm")
# positive relationship but its not crazy

# plot of urbanization and number of checklists
ggplot(summary_filt, aes(x=urban, y=log(number_checklists)))+
  geom_point()+
  geom_smooth(method="lm")
# need to check this outlier with 50,000 checklists

# plot of urbanization and latitude
ggplot(summary_filt, aes(x=y, y=urban))+
  geom_point()+
  geom_smooth()
# no real pattern

# plot of SR and latitude
ggplot(summary_filt, aes(x=y, y=log(total_SR)))+
  geom_point()+
  geom_smooth() +
  facet_wrap(~urban)
# peaks at intermediate latitudes (good)

ggplot(summary_filt, aes(x=abs(y), y=log(total_SR)))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~urban)

ggplot(summary_filt, aes(x=x, y=log(total_SR)))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~urban)

ggplot(summary_filt, aes(x=x, y=y))+
  geom_point()+
  facet_wrap(~urban)

#########################################
# Re-do extraction of the urbanization values with the raster that I removed values from

# load raster with the high impact low urbanization areas filtered out
GHSL_filt <- rast("SMOD_global/urbanization.tif")

# load csv of filtered summary data
GHSL_filt$GHS_SMOD_E2020_GLOBE_R2022A_54009_1000_V1_0
# extract the values from this raster into another column

#summary_filt$urban_filtered <- as.data.frame(terra::extract(GHSL_filt, summary_filt[,c(20:21)]))$GHS_SMOD_E2020_GLOBE_R2022A_54009_1000_V1_0

#summary_new_filt <- summary_filt %>% filter(!is.na(urban_filtered))
# that removed almost another 10,000

#write.csv(summary_new_filt, "5yr_summary/summary_thresholded.csv")

summary_new_filt <- read.csv("5yr_summary/summary_thresholded.csv")

# plot of urbanization and species richness
ggplot(summary_new_filt, aes(x=urban, y=total_SR))+
  geom_point()+
  geom_smooth(method="lm")
# definitely going down which is great

# plot of species richness and number of checklists
ggplot(summary_new_filt, aes(x=log(number_checklists), y=total_SR))+
  geom_point()+
  geom_smooth()
# positive relationship but its not crazy

summary_new_filt %>% 
  filter(number_checklists < 10000) %>% # filter for checklists less than 10000
  ggplot(aes(x=number_checklists, y=total_SR))+
  geom_point()+
  geom_smooth()

# plot of urbanization and number of checklists
ggplot(summary_new_filt, aes(x=urban, y=log(number_checklists)))+
  geom_point()+
  geom_smooth(method="lm")
# need to check this outlier with 50,000 checklists


ggplot(summary_new_filt, aes(x=abs(y), y=log(number_checklists)))+
  geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~urban_filtered)


# plot of urbanization and latitude
ggplot(summary_new_filt, aes(x=y, y=urban))+
  geom_point()+
  geom_smooth()
# no real pattern

# plot of SR and latitude
ggplot(summary_filt, aes(x=y, y=total_SR))+
  geom_point()+
  geom_smooth()
# peaks at intermediate latitudes (good)

# plot of absolute latitude and species richness
ggplot(summary_new_filt, aes(x=abs(y), y=total_SR))+
  geom_point()+
  geom_smooth()


ggplot(summary_new_filt, aes(x=abs(y), y=total_SR))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~urban)

ggplot(summary_new_filt, aes(x=x, y=total_SR))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~urban)

ggplot(summary_new_filt, aes(x=x, y=y))+
  geom_point()+
  facet_wrap(~urban)

hist(summary_new_filt$total_SR)
# does not need to be logged



####################
# Apply higher threshold (95th quartile)

summary_new_filt %>% filter(number_checklists>=95) %>% 
  ggplot(aes(x=abs(y), y=total_SR))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~urban)
# definitely looks less strong in the urban areas


###########################
# Work on figure with urbanization and latitudinal gradient
summary_new_filt <- read.csv("5yr_summary/summary_thresholded.csv")
# want to tun into 3 levels, urban, peri-urban, and urban center
urb_levels <- summary_new_filt

urb_levels$urban[which(urb_levels$urban %in% c(11, 12, 13))] <- 1
urb_levels$urban[which(urb_levels$urban %in% c(21, 22, 23))] <- 2
urb_levels$urban[which(urb_levels$urban==30)] <- 3

urban_names <- c(
  `1` = "Rural",
  `2` = "Peri-Urban",
  `3` = "Urban"
)

LDG <- ggplot(urb_levels, aes(x=abs(y), y=total_SR))+
  geom_point(alpha=0.2, shape=1)+
  geom_smooth(method="lm", color="forestgreen")+
  labs(x="Absolute Latitude", y="Species Richness")+
  facet_wrap(~urban, labeller = labeller(urban = urban_names))+
  theme_bw()

# save as pdf
pdf("urban_LDG.pdf", height=5, width=7)
LDG
dev.off()

# try with overall latitude (and a quadratic model)
ggplot(urb_levels, aes(x=y, y=total_SR))+
  geom_point(alpha=0.2, shape=1)+
  geom_smooth(method="lm", formula = y ~ x + I(x^2), color="forestgreen")+
  labs(x="Absolute Latitude", y="Species Richness")+
  facet_wrap(~urban, labeller = labeller(urban = urban_names))+
  theme_bw()



##############################################
# Rasterize and plot data
# load species richness data
summary_filt <- read.csv("5yr_summary/summary_thresholded.csv")
# load urbanization raster
GHSL_filt <- rast("SMOD_global/urbanization.tif")

#summary_vect <- vect(summary_filt, crs=crs(GHSL_filt), geom=c("x","y"))
#thresholded.raster <- terra::rasterize(summary_vect, GHSL_filt, field="total_SR") 
#writeRaster(thresholded.raster, "5yr_summary/summary_rasterized.tif", overwrite=TRUE)
# plot by each bbox

thresholded.raster <- rast("5yr_summary/summary_rasterized.tif")
plot(thresholded.raster)
# load bbox extents
bbox <- read.csv("bounding_box_coordinates.csv")
bbox <- bbox[-c(7:9),] # remove the cut up r2c2
bbox[6,2]<-"r2c2"
bbox[6,4]<- 0 # replacing xmin




names <- bbox$names

##### Plot with points and species richness
for (i in 1:nrow(bbox)){
  xlim <- c(bbox$xmin[i], bbox$xmax[i])
  ylim <- c(bbox$ymin[i], bbox$ymax[i])
  extent <- c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])
  #cropped <- crop(thresholded.raster, ext(extent))
  
  # try with points to be able to actually see things
  assign(paste0("plot_", bbox$names[i]),
         ggplot(data=world)+
           geom_sf() +
           geom_point(data=summary_filt, aes(x=x, y=y, color=total_SR), size=0.1) +
           coord_sf(crs = 4326, xlim = xlim, ylim = ylim, expand = FALSE) +
           labs(x="Longitude", y="Latitude", color="Richness")+
           scale_color_viridis_c(na.value = NA)+
           theme_bw())
}
# save plots as pdf
# make plot of whole globe
world_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_filt, aes(x=x, y=y, color=total_SR), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  labs(x="Longitude", y="Latitude", color="Richness")+
  theme_bw()

pdf("point_plots.pdf")
world_plot
plot_r1c1
plot_r1c2
plot_r1c3
plot_r1c4
plot_r2c1
plot_r2c2
plot_r2c3
plot_r2c4
plot_r3c1
plot_r3c2
plot_r3c4
plot_r4c1
plot_r4c2
plot_r4c3
plot_r4c4
dev.off()


######### Plots with the raster data

for (i in 1:nrow(bbox)){
  xlim <- c(bbox$xmin[i], bbox$xmax[i])
  ylim <- c(bbox$ymin[i], bbox$ymax[i])
  extent <- c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])
  cropped <- crop(thresholded.raster, ext(extent))
  
  # try with points to be able to actually see things
  assign(paste0("rast_", bbox$names[i]),
         ggplot(data=world)+
           geom_sf() +
           stat_spatraster(data=cropped, aes()) +
           coord_sf(crs = 4326, xlim = xlim, ylim = ylim, expand = FALSE) +
           labs(x="Longitude", y="Latitude", color="Richness") +  
           scale_fill_viridis_c(na.value = NA)+
           theme_bw())
}
# save plots as pdf
# make plot of whole globe
world_rast<-ggplot(data=world)+
  geom_sf() +
  stat_spatraster(data=thresholded.raster, aes()) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_fill_viridis_c(na.value = NA)+
  theme_bw()

pdf("rast_plots.pdf")
world_rast
rast_r1c1
rast_r1c2
rast_r1c3
rast_r1c4
rast_r2c1
rast_r2c2
rast_r2c3
rast_r2c4
rast_r3c1
rast_r3c2
rast_r3c4
rast_r4c1
rast_r4c2
rast_r4c3
rast_r4c4
dev.off()
# these do not look very good


########## Plots with points colored by number of checklists
# remove checklist outliers
summary_new <- summary_filt %>% filter(number_checklists<10000)
for (i in 1:nrow(bbox)){
  xlim <- c(bbox$xmin[i], bbox$xmax[i])
  ylim <- c(bbox$ymin[i], bbox$ymax[i])
  extent <- c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])
  #cropped <- crop(thresholded.raster, ext(extent))
  
  # try with points to be able to actually see things
  assign(paste0("check_", bbox$names[i]),
         ggplot(data=world)+
           geom_sf() +
           geom_point(data=summary_new, aes(x=x, y=y, color=number_checklists), size=0.1) +
           coord_sf(crs = 4326, xlim = xlim, ylim = ylim, expand = FALSE) +
           scale_color_viridis_c(na.value = NA, option="C")+
           labs(x="Longitude", y="Latitude", color="Checklists")+
           theme_bw())
}
# save plots as pdf
# make plot of whole globe
world_check <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_new, aes(x=x, y=y, color=number_checklists), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="D")+
  labs(x="Longitude", y="Latitude", color="Checklists")+
  theme_bw()

pdf("checklist_plots.pdf")
world_check
check_r1c1
check_r1c2
check_r1c3
check_r1c4
check_r2c1
check_r2c2
check_r2c3
check_r2c4
check_r3c1
check_r3c2
check_r3c4
check_r4c1
check_r4c2
check_r4c3
check_r4c4
dev.off()


######## Plots with urbanization scores
for (i in 1:nrow(bbox)){
  xlim <- c(bbox$xmin[i], bbox$xmax[i])
  ylim <- c(bbox$ymin[i], bbox$ymax[i])
  extent <- c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])
  #cropped <- crop(thresholded.raster, ext(extent))
  
  # try with points to be able to actually see things
  assign(paste0("urban_", bbox$names[i]),
         ggplot(data=world)+
           geom_sf() +
           geom_point(data=summary_filt, aes(x=x, y=y, color=urban), size=0.1) +
           coord_sf(crs = 4326, xlim = xlim, ylim = ylim, expand = FALSE) +
           labs(x="Longitude", y="Latitude", color="Urban")+
           scale_color_viridis_c(na.value = NA, option="F")+
           theme_bw())
}
# save plots as pdf
# make plot of whole globe
world_urban <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_filt, aes(x=x, y=y, color=urban), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="F")+
  labs(x="Longitude", y="Latitude", color="Urban")+
  theme_bw()

pdf("urban_plots.pdf")
world_urban
urban_r1c1
urban_r1c2
urban_r1c3
urban_r1c4
urban_r2c1
urban_r2c2
urban_r2c3
urban_r2c4
urban_r3c1
urban_r3c2
urban_r3c4
urban_r4c1
urban_r4c2
urban_r4c3
urban_r4c4
dev.off()


####### Put four plots together for writing class figures
library(ggpubr)
plot_arrange <- ggarrange(world_plot, world_check, labels=c("A", "B"), 
                          nrow=2, ncol=1, label.x=0.1)

ggsave(plot=plot_arrange, filename="point_maps.pdf")





###### Plots of species richness but only in highly urbanized areas
summary_urb <- summary_filt %>% filter(urban==30)
for (i in 1:nrow(bbox)){
  xlim <- c(bbox$xmin[i], bbox$xmax[i])
  ylim <- c(bbox$ymin[i], bbox$ymax[i])
  extent <- c(bbox$xmin[i], bbox$xmax[i], bbox$ymin[i], bbox$ymax[i])
  #cropped <- crop(thresholded.raster, ext(extent))
  
  # try with points to be able to actually see things
  assign(paste0("only_urb_", bbox$names[i]),
         ggplot(data=world)+
           geom_sf() +
           geom_point(data=summary_urb, aes(x=x, y=y, color=total_SR), size=0.1) +
           coord_sf(crs = 4326, xlim = xlim, ylim = ylim, expand = FALSE) +
           labs(x="Longitude", y="Latitude", color="Richness")+
           scale_color_viridis_c(na.value = NA)+
           theme_bw())
}
# save plots as pdf
# make plot of whole globe
world_only_urb <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_urb, aes(x=x, y=y, color=total_SR), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  labs(x="Longitude", y="Latitude", color="Richness")+
  theme_bw()

pdf("only_urb_plots.pdf")
world_only_urb
only_urb_r1c1
only_urb_r1c2
only_urb_r1c3
only_urb_r1c4
only_urb_r2c1
only_urb_r2c2
only_urb_r2c3
only_urb_r2c4
only_urb_r3c1
only_urb_r3c2
only_urb_r3c4
only_urb_r4c1
only_urb_r4c2
only_urb_r4c3
only_urb_r4c4
dev.off()




# points colored by number of checklists
ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_filt, aes(x=x, y=y, color=log(number_checklists)), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  theme_bw()


# points colored by urbanization score
ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_filt, aes(x=x, y=y, color=urban), size=0.1, alpha=0.5) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  theme_bw()
# some good data from highly urbanized areas


# only highly urbanized
high_urban <- summary_filt %>% filter(urban==30)
ggplot(data=world)+
  geom_sf() +
  geom_point(data=high_urban, mapping=aes(x=x, y=y), size=0.1, alpha=0.5, color="red") +
  coord_sf(crs = 4326, expand = FALSE) +
  # scale_color_viridis_c(na.value = NA)+
  theme_bw()


# only low urbanization
low_urban <- summary_filt %>% filter(urban%in% c(11, 12, 13))
ggplot(data=world)+
  geom_sf() +
  geom_point(data=low_urban, mapping=aes(x=x, y=y), size=0.1, alpha=0.5, color="forestgreen") +
  coord_sf(crs = 4326, expand = FALSE) +
  # scale_color_viridis_c(na.value = NA)+
  theme_bw()


# summary of rasters by urbanization level
summary_filt %>% group_by(urban) %>% summarise(n=n())
# definitely unbalanced but not too bad actually

ggplot(data=world)+
  geom_sf() +
  geom_point(data=low_urban, mapping=aes(x=x, y=y), size=0.1, alpha=0.5, color="forestgreen") +
  coord_sf(crs = 4326, expand = FALSE) +
  # scale_color_viridis_c(na.value = NA)+
  theme_bw()+
  facet_wrap(~urban)


#####################################################
# Try to figure out what cities the urban sites are in

# load shapefile of cities
cities <- read.csv("World_Cities.csv")
# turn into spatvector
cities_vect <- vect(cities, crs=crs(GHSL_filt), geom=c("X","Y"))
# turn data with high urbanization into spatvector
urb_vect <- vect(summary_urb, crs=crs(GHSL_filt), geom=c("x","y"))

urb_cities <- as.data.frame(nearest(urb_vect, cities_vect))

# merge to get city names
cities_in_df <-merge(urb_cities, cities, by.x="to_id", by.y="FID", all.x=TRUE)

length(unique(cities_in_df$CITY_NAME))
# 621 cities
# not all of these are totally accurate because they don't have some smaller cities 
 # (e.g. labelled bellingham as vancouver)
# but gives somewhat of an idea

length(unique(cities_in_df$CNTRY_NAME))
# 121 countries

countries <- cities_in_df %>% group_by(CNTRY_NAME) %>% summarise(n=n())
cities <- cities_in_df %>% group_by(CITY_NAME) %>% summarise(n=n())


##############################
# Classifying points into biomes using a terrestrial biomes shapefile
biomes <- st_read("/Volumes/Expansion/eBird/wwf_biomes/wwf_terr_ecos.shp")

biomes_vect <- vect(biomes)

summary_vect <- vect(summary_filt, crs=crs(biomes_vect), geom=c("x","y"))
  
biomes_overlap_df <- as.data.frame(intersect(summary_vect, biomes_vect))

biomes_overlap[[1]]
head(biomes_overlap_df)

#############################
# Coverage of points in urban, peri-urban, and rural areas

summary_filt <- read.csv("5yr_summary/summary_thresholded.csv")
urb_levels <- summary_filt

urb_levels$urban[which(urb_levels$urban %in% c(11, 12, 13))] <- "Rural"
urb_levels$urban[which(urb_levels$urban %in% c(21, 22, 23))] <- "Peri-Urban"
urb_levels$urban[which(urb_levels$urban==30)] <- "Urban"

### Thresholded at 53 (mean of checklists needed for 95% coverage)
threshold_50 <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=urb_levels, aes(x=x, y=y, color=total_SR), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  facet_wrap(~urban, nrow=2)+
  theme_bw() +
  labs(x="Longitude", y="Latitude", color="Richness", title="Thresholded at mean (53 checklists)",
       caption = "n=75988")

ggsave(filename="threshold_50_coverage.png", plot=threshold_50, height=8, width=10)

# plot of LDG
LDG_50 <- ggplot(urb_levels, aes(x=abs(y), y=total_SR))+
  geom_point(alpha=0.2, shape=1)+
  geom_smooth(method="lm", color="forestgreen")+
  labs(x="Absolute Latitude", y="Species Richness", title="Thresholded at 53 checklists")+
  facet_wrap(~urban)+
  theme_bw()


##### Thresholded at 67 (75% quantile)
urb_levels_75 <- urb_levels %>% filter(number_checklists>=67) # 62466
threshold_75 <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=urb_levels_75, aes(x=x, y=y, color=total_SR), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  facet_wrap(~urban)+
  labs(x="Longitude", y="Latitude", color="Richness", title="Thresholded at 75 (67 checklists)",
       caption = "n=62466") +
  theme_bw()
ggsave(filename="threshold_75_coverage.png", plot=threshold_75, height=8, width=10)

# Plot of LDG
LDG_75 <- ggplot(urb_levels_75, aes(x=abs(y), y=total_SR))+
  geom_point(alpha=0.2, shape=1)+
  geom_smooth(method="lm", color="forestgreen")+
  labs(x="Absolute Latitude", y="Species Richness", title="Thresholded at 67 checklists")+
  facet_wrap(~urban)+
  theme_bw()


##### Thresholded at 95 (95th percentile)
urb_levels_95 <- urb_levels %>% filter(number_checklists>=95) # 62466
threshold_95 <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=urb_levels_95, aes(x=x, y=y, color=total_SR), size=0.1) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  facet_wrap(~urban)+
  labs(x="Longitude", y="Latitude", color="Richness", title="Thresholded at 95 (95 checklists)",
       caption = "n=46161") +
  theme_bw()
ggsave(filename="threshold_95_coverage.png", plot=threshold_95, height=8, width=10)

# LDG
LDG_95 <- ggplot(urb_levels_75, aes(x=abs(y), y=total_SR))+
  geom_point(alpha=0.2, shape=1)+
  geom_smooth(method="lm", color="forestgreen")+
  labs(x="Absolute Latitude", y="Species Richness", title="Thresholded at 95 checklists")+
  facet_wrap(~urban)+
  theme_bw()

# they all look pretty similar

# put LDG plots together
ggarrange(LDG_50, LDG_75, LDG_95)


