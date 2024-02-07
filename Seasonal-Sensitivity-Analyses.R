########## Sensitivity analyses for seasonal data #########
library(tidyverse)
library(terra)
library(marginaleffects)
library(sf)
library(elevatr)


##### 1) Larger and smaller cutoff for number of checklists
## Threshold of 90% coverage for winter is 38 checklists and for summer is 33 checklists
## Threshold of 98% coverage winter is 175 checklists and for summer is 144 checklists

## I'll start with higher coverage (98%) because I can just filter it from the modeling data
model.dat <- read.csv("season_modeling_data.csv")
wint <- model.dat %>% filter(season=="Winter" & number_checklists >= 175)
sum <- model.dat %>% filter(season=="Summer" & number_checklists >= 144)
dat.98 <- rbind(wint, sum)

# model 
season.model.98 <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                        precip + log(number_checklists) + elevation, dat.98)

square <- function(x){
  x^2
}
predicted.season.98 <- avg_predictions(season.model.98, by=c("abslat", "urban2", "season"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                       season = c("Summer", "Winter")))


#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot.98 <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat.98, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season.98, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season.98, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot.98
# yup pretty much the same result
plot_slopes(season.model.98, variables="abslat", condition=c("urban2", "season"))
# yes, the slope overlaps 0 in summer urban still


###### Now run lower coverage 
summer_dat <- read.csv("summer_richness_summary.csv") # full summary dataset
summer_dat <- summer_dat %>% filter(number_checklists >= 33)
summer_dat$season <- "Summer"
winter_dat <- read.csv("winter_richness_summary.csv")
winter_dat <- winter_dat %>% filter(number_checklists >= 38)
winter_dat$season <- "Winter"
season_dat <- rbind(summer_dat, winter_dat)

# now add in other stuff

# extract urbanization values
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtMollweide.tif")
season_dat$x <- xFromCell(GHSL, season_dat$cell) # extract the coordinates from the cells
season_dat$y <- yFromCell(GHSL, season_dat$cell)
season_dat$urban <- as.data.frame(terra::extract(GHSL, season_dat[,c(10:11)]))$SMOD_global
season_dat.90 <- season_dat %>% filter(!is.na(urban))
# Extract elevation
summary_sf <- st_as_sf(season_dat.90, coords=c("x", "y"), crs=st_crs(GHSL))
dat_latlong <- st_transform(summary_sf, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                     lat = sf::st_coordinates(.)[,2])
latlong_df <- get_elev_point(latlong_df[,c(12,13,2:10)], prj=crs(summary_sf), src="aws", overwrite=TRUE) # extract elevations from amazon web services
season_dat.90 <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                             y = sf::st_coordinates(.)[,2]))
# add precipitation
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
season_dat.90$precip <- as.data.frame(terra::extract(precip, season_dat.90[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12
# add hemisphere
season_dat.90$hemisphere <- "northern"
season_dat.90$hemisphere[season_dat.90$lat<0]<-"southern"
# add abslat
season_dat.90$abslat <- abs(season_dat.90$lat)
## turn urban into 3 categories
season_dat.90 <- season_dat.90 %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
season_dat.90$urban2 <- as.factor(season_dat.90$urban2) # turn into factor for modelling
season_dat.90$urban2 <- factor(season_dat.90$urban2, levels = c("1", "2", "3"),
                               labels = c("Natural", "Suburban", "Urban")) # rename as natural, suburban, and urban


### Model
season.model.90 <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                        precip + log(number_checklists) + elevation, season_dat.90)

square <- function(x){
  x^2
}
predicted.season.90 <- avg_predictions(season.model.90, by=c("abslat", "urban2", "season"), transform=square, 
                                       newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                          season = c("Summer", "Winter")))


#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot.90 <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(season_dat.90, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season.90, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season.90, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot.90
# yup pretty much the same result
plot_slopes(season.model.90, variables="abslat", condition=c("urban2", "season"))
# yes, the slope overlaps 0 in summer urban still!!








##### 2) Different thresholds of modification
lowmod <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtLowThreshold.tif")
highmod <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtHighThreshold.tif")



## Do highmod first because I can just filter it from the modeling data (this is removing everything with more than values of 0.25 HI)
dat_highmod <- read.csv("season_modeling_data.csv")

# extract raster data
dat_highmod$urban <- as.data.frame(terra::extract(highmod, dat_highmod[,c(18:19)]))$SMOD_global
# remove NA values
dat_highmod <- dat_highmod %>% filter(!is.na(urban))

### Model
season.highmod <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                        precip + log(number_checklists) + elevation, dat_highmod)

square <- function(x){
  x^2
}
predicted.season.highmod <- avg_predictions(season.highmod, by=c("abslat", "urban2", "season"), transform=square, 
                                       newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                          season = c("Summer", "Winter")))


#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot.highmod <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat_highmod, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season.highmod, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season.highmod, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot.highmod
# yup pretty much the same result
plot_slopes(season.highmod, variables="abslat", condition=c("urban2", "season"))
# yes, the slope overlaps 0 in summer urban still!!



#### Now with lowmod
summer_dat <- read.csv("summer_richness_summary.csv")
summer_filt95 <- summer_dat %>% filter(number_checklists >= 62) 
winter_dat <- read.csv("winter_richness_summary.csv")
winter_filt95 <- winter_dat %>% filter(number_checklists >= 73) 
season_lowmod <- rbind(summer_filt95, winter_filt95)


season_dat$x <- xFromCell(GHSL, season_dat$cell) # extract the coordinates from the cells
season_dat$y <- yFromCell(GHSL, season_dat$cell)
season_dat$urban <- as.data.frame(terra::extract(lowmod, season_dat[,c(10:11)]))$SMOD_global
season_dat.lowmod <- season_dat %>% filter(!is.na(urban))
# Extract elevation
summary_sf <- st_as_sf(season_dat.lowmod, coords=c("x", "y"), crs=st_crs(lowmod))
dat_latlong <- st_transform(summary_sf, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                     lat = sf::st_coordinates(.)[,2])
latlong_df <- get_elev_point(latlong_df[,c(12,13,2:10)], prj=crs(summary_sf), src="aws", overwrite=TRUE) # extract elevations from amazon web services
season_dat.lowmod <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                                  y = sf::st_coordinates(.)[,2]))
# add precipitation
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
season_dat.lowmod$precip <- as.data.frame(terra::extract(precip, season_dat.lowmod[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12
# add hemisphere
season_dat.lowmod$hemisphere <- "northern"
season_dat.lowmod$hemisphere[season_dat.lowmod$lat<0]<-"southern"
# add abslat
season_dat.lowmod$abslat <- abs(season_dat.lowmod$lat)
## turn urban into 3 categories
season_dat.lowmod <- season_dat.lowmod %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
season_dat.lowmod$urban2 <- as.factor(season_dat.lowmod$urban2) # turn into factor for modelling
season_dat.lowmod$urban2 <- factor(season_dat.lowmod$urban2, levels = c("1", "2", "3"),
                               labels = c("Natural", "Suburban", "Urban")) # rename as natural, suburban, and urban

## Model

### Model
season.lowmod <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                       precip + log(number_checklists) + elevation, season_dat.lowmod)

square <- function(x){
  x^2
}
predicted.season.lowmod <- avg_predictions(season.highmod, by=c("abslat", "urban2", "season"), transform=square, 
                                            newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                               season = c("Summer", "Winter")))


#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot.lowmod <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(season_dat.lowmod, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season.lowmod, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season.lowmod, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot.lowmod
# yup pretty much the same result
plot_slopes(season.lowmod, variables="abslat", condition=c("urban2", "season"))
# yes, the slope overlaps 0 in summer urban still!!




