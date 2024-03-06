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


# try to rasterize to plot
# try to aggregate to make more visible
#plot(thresholded.raster)
#raster.aggregated <- aggregate(thresholded.raster, fact = 10, fun="sum")
#raster.aggregated[is.na(raster.aggregated)] <- 10000
#plot(raster.aggregated)



# load summary data
summary <- read.csv("global_richness_summary.csv")
# 2075417 squares with richness values (this is the unfiltered data)

# add the lat and long values

GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
plot(GHSL)

summary$x <- xFromCell(GHSL, summary$cell) # extract the coordinates from the cells
summary$y <- yFromCell(GHSL, summary$cell)

summary$urban <- as.data.frame(terra::extract(GHSL, summary[,c(21:22)]))$SMOD_global

# plot data
summary_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA)+
  labs(x="Longitude", y="Latitude", color="SR", title="All data included")+
  theme_bw()

# this is the plot with nothing filtered out, very good coverage




###### Load and plot filtered data
summary_filt <- read.csv("5yr_summary/summary_thresholded.csv")

summary95_plot <- ggplot(data=world)+
  geom_sf(lwd=2) +
  geom_point(data=summary_filt, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="C")+
  labs(x="Longitude", y="Latitude", color="Species \nRichness")+
  theme_bw()
ggsave(summary95_plot, file="CoveragePlot.png", height=8, width=6)

####### Try the other filters and plot to compare
# 97
summary_filt97 <- summary_filt %>% filter(number_checklists >= 139)
summary97_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_filt97, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude", color="SR", title="95 of 97% coverage")+
  theme_bw()

# 98
summary_filt98 <- summary_filt %>% filter(number_checklists >= 203)
summary98_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summary_filt98, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude", color="SR", title="95 of 98% coverage")+
  theme_bw()

# save as a png
pdf("coverage_plots.pdf")
summary_plot
summary95_plot
summary97_plot
summary98_plot
dev.off()


####################
# Some plots of relationships between variables

# plot of urbanization and species richness
ggplot(summary_filt, aes(x=urban, y=total_SR))+
  geom_point()+
  geom_smooth(method="lm")
# species richness decreases with urbanization, good!

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


###########################
# Plot latitudinal gradient and urbanization

# want to tun into 3 levels, urban, peri-urban, and urban center
urb_levels <- summary_filt

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
  theme_bw() # has the expected relationship!

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
# Plot data as points in different squares

# load bbox extents
bbox <- read.csv("bounding_box_coordinates.csv")
bbox <- bbox[-c(7:9),] # remove the cut up r2c2
bbox[6,2]<-"r2c2"
bbox[6,4]<- 0 # replacing xmin

names <- bbox$names


##### Plots with points and species richness in different areas
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

pdf("point_plots.pdf")
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



# summary of rasters by urbanization level
summary_filt %>% group_by(urban) %>% summarise(n=n())
# definitely unbalanced but not too bad actually



#####################################################
# LOOK AT WHAT CITIES URBAN AREAS ARE IN
dat.urb <- read.csv("modeling_data.csv") %>% filter(urban2=="Urban")
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif")
# load cities - just point data
cities <- read.csv("/Volumes/Expansion/eBird/World_Cities.csv")
# turn into spatvector
cities_vect <- vect(cities, crs=crs(GHSL), geom=c("X","Y"))
# turn data with high urbanization into spatvector
urb_vect <- vect(dat.urb, crs=crs(GHSL), geom=c("long","lat"))

urb_cities <- as.data.frame(nearest(urb_vect, cities_vect))

# merge to get city names
cities_in_df <-merge(urb_cities, cities, by.x="to_id", by.y="FID", all.x=TRUE)

length(unique(cities_in_df$CITY_NAME))
# 624 cities
# not all of these are totally accurate because they don't have some smaller cities 
 # (e.g. labelled bellingham as vancouver)
# but gives somewhat of an idea

length(unique(cities_in_df$CNTRY_NAME))
# 123 countries

countries <- cities_in_df %>% group_by(CNTRY_NAME) %>% summarise(n=n())
cities <- cities_in_df %>% group_by(CITY_NAME) %>% summarise(n=n())




### Maps with urbanization score
dat <- read.csv("modeling_data.csv")
world <- ne_countries(scale = "medium", returnclass = "sf")

urb.map<-ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=dat, aes(x=long, y=lat, color=urban2), size=0.05, alpha=0.3) +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = FALSE) +
  labs(x="Longitude", y="Latitude")+
 # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
#  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_void()+
  theme(legend.title=element_blank(), legend.position = c(.85, .3), text=element_text(size=15), strip.text=element_blank())+
  facet_wrap(~urban2, ncol=2)
urb.map
# there are like 2 points in the arctic zone
ggsave(urb.map, file="urban.map.png", height=5, width=10)

# try all together
map.together <- ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=dat, aes(x=long, y=lat, color=urban2), size=0.05, alpha=0.3) +
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = FALSE) +
  labs(x="Longitude", y="Latitude")+
 # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
#  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_void() +
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=15))
ggsave(map.together, file="map.together.png", height=5, width=7)


geographic_zones <- 
  ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
#  geom_point(data=dat, aes(x=long, y=lat, color=urban2), size=0.05) +
 # scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  coord_sf(crs = 4326, expand = FALSE, ylim=c(-65,70)) +
  labs(x="Longitude", y="Latitude")+
  geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50), alpha=0.7, lty=3)+ # geographic zones
  geom_hline(yintercept=0, alpha=0.8, lty=2) + # for equator
  theme_void()+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-50, fill="lightblue", alpha=0.2)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=50, ymax=Inf, fill="lightblue", alpha=0.2)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=35, ymax=50, fill="forestgreen", alpha=0.2)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-50, ymax=-35, fill="forestgreen", alpha=0.2)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=23.4, ymax=35, fill="yellow", alpha=0.2)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-35, ymax=-23.4, fill="yellow", alpha=0.2)+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-23.4, ymax=23.4, fill="orange", alpha=0.2)+
  theme(legend.title=element_blank(), legend.position = c(.85, .15), text=element_text(size=15))
ggsave(geographic_zones,  file="zones_map.png")


#library(maps)
#map('world',col="white", fill=TRUE, bg="white", lwd=0.2, mar=rep(0,4), border=1, ylim=c(-80,80))  




