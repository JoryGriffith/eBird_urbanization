###### Writing code for sensitivity analyses ##########
library(tidyverse)
library(terra)
library(marginaleffects)
library(sf)
library(elevatr)
library(ggeffects)

##### 1) Larger and smaller cutoff for number of checklists
## Threshold of 90% coverage is 44 checklists
## Threshold of 98% coverage is 198 checklists


## I'll start with higher coverage because I can just filter it from the modeling data
dat.98 <- read.csv("modeling_data.csv") %>% filter(number_checklists >= 198) # 30,736
dat.98 %>% filter(CONTINENT=="Antarctica") # there is one point in antarctica, latitude is -54.05 so it just misses the cutoff
dat.98 <- dat.98 %>% filter(lat <= 70 & lat >=-55)
# neither of the filters did anything
model.98 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat.98)
summary(model.98)
square <- function(x){
  x^2
} 

## plot results
predicted.98<-avg_predictions(model.98, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))


LDG98.plot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat.98, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.98, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.98, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
LDG98.plot
# yup looks pretty much the same!
plot_slopes(model.98, variables="abslat", condition=c("urban2"))


######## Now thinned
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL)) # sample grid

# assign cell number to each point in my data
vect <- st_as_sf(dat.98, crs=st_crs(GHSL), coords=c("x","y"))
xy=st_coordinates(vect)
# get cell number that each point is in
dat.98$cell.subsample<-cellFromXY(sample.grid, xy)


square <- function(x){
  x^2
} # make function to square
predicted <- list()
ggeffects.slopes <- list()
ggeffects.slopes.contrast <- list()
set.seed(30)

for (i in 1:1000){
  dat.thinned <- dat.98 %>% group_by(cell.subsample, urban2) %>% sample_n(1) # randomly sample point from each vell
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned) # run model
  predicted[[i]] <- avg_predictions(lm.thinned, by=c("abslat", "urban2"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                       urban2=c("Natural", "Suburban", "Urban"))) # store predictions for plotting
  ggeffects.slopes[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"), test = NULL) # store slopes
    ggeffects.slopes.contrast[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2")) # store contrasts
}

predicted_df <- bind_rows(predicted)
write.csv(predicted_df, "supplement_figs/98.coverage.thinned.results.csv")

ggeffects.slopes.df <- bind_rows(ggeffects.slopes)
write.csv(ggeffects.slopes.df, file="supplement_figs/98.coverage.thinned.slopes.csv")
ggslopes <- ggeffects.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), conf.low=min(conf.low))
ggslopes


ggeffects.contrast_df <- bind_rows(ggeffects.slopes.contrast)
#write.csv(ggeffects.contrast_df, file="thinned_results/thinned_contrasts.csv")
contrast <- ggeffects.contrast_df %>% filter(urban2=="Suburban-Urban") %>% filter(p.value<0.05) # 989
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Suburban") %>% filter(p.value<0.05) #1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Urban") %>% filter(p.value<0.05) # 1000

######## Plot results
thinned.results.98 <- read.csv("supplement_figs/98.coverage.thinned.results.csv")
thinned.results.summary.98 <- thinned.results.98 %>% group_by(abslat, urban2) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))

thinned.plots.98 <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat.98, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary.98, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(thinned.results.summary.98, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,600))+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=18), axis.title=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
# this is the plot with the 95% of the confidence intervals
thinned.plots.98
#ggsave(thinned.plots.98, file="supplement_figs/thinned.results.coverage98.png", height=6, width=7)



##############################
###############################




##### Now with threshold of 90% coverage ################
summary <- read.csv("global_richness_summary.csv") # load full summary dataset
summary_filt90 <- summary %>% filter(number_checklists >= 44) # 143,762
# merge with urbanization raster values
GHSL<-rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtMollweide.tif")

# extract cell numbers, urbanization scores, and x and y coordinates from the raster
summary_filt90$x <- xFromCell(GHSL, summary_filt90$cell) # extract the coordinates from the cells
summary_filt90$y <- yFromCell(GHSL, summary_filt90$cell)

summary_filt90$urban <- as.data.frame(terra::extract(GHSL, summary_filt90[,c(21:22)]))$SMOD_global

summary_filt90 %>% group_by(urban) %>% summarise(n=n())

# remove ones with NaN urbanization score
summary_filt90 <- summary_filt90 %>% filter(!is.na(urban)) # 110,211 datapoints now


# Turn into sf (spatial) object
summary_sf <- st_as_sf(summary_filt90, coords=c("x", "y"), crs=st_crs(GHSL))

## Add elevation
# need to convert the points back to lat long
summary_latlong <- st_transform(summary_sf, crs=st_crs(4326))
latlong_df <- summary_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                     lat = sf::st_coordinates(.)[,2])

crs(summary_latlong)
latlong_df <- get_elev_point(latlong_df[,c(23,24,2:21)], prj=crs(summary_latlong), src="aws", overwrite=TRUE) # extract elevations from amazon web services
# does the crs change anything
#test <- latlong_df
#test$elevation <- as.data.frame(get_elev_point(test, prj=crs(GHSL), src="aws", overwrite=TRUE))[,1]

dat_summary90 <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                             y = sf::st_coordinates(.)[,2]))
dat_summary90

#### assign hemisphere
dat_summary90$hemisphere <- "northern"
dat_summary90$hemisphere[dat_summary90$lat<0]<-"southern"


## add precipitation
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
dat_summary90$precip <- as.data.frame(terra::extract(precip, dat_summary90[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12

dat_summary90$abslat <- abs(dat_summary90$lat)

## turn urban into 3 categories
dat_summary90 <- dat_summary90 %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat_summary90 %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat_summary90$urban2 <- as.factor(dat_summary90$urban2) # turn into factor for modelling

dat_summary90$urban2 <- factor(dat_summary90$urban2, levels = c("1", "2", "3"),
                     labels = c("Natural", "Suburban", "Urban")) # rename as natural, suburban, and urban





###### Model
model90 <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + log(number_checklists) +
                   precip + elevation, dat_summary90)
anova(model90)

square <- function(x){
  x^2
} 

## plot results
predicted.full<-avg_predictions(model90, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                   urban2=c("Natural", "Suburban", "Urban")))


LDG90.plot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat_summary90, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
LDG90.plot

# these both change the lines a bit but the takeaways are still the same
plot_slopes(model90, variables="abslat", condition=c("urban2"))


########### THIN
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL)) # sample grid

# assign cell number to each point in my data
vect <- st_as_sf(dat_summary90, crs=st_crs(GHSL), coords=c("x","y"))
xy=st_coordinates(vect)
# get cell number that each point is in
dat_summary90$cell.subsample<-cellFromXY(sample.grid, xy)


square <- function(x){
  x^2
} # make function to square
predicted <- list()
ggeffects.slopes <- list()
ggeffects.slopes.contrast <- list()
set.seed(30)

for (i in 1:1000){
  dat.thinned <- dat_summary90 %>% group_by(cell.subsample, urban2) %>% sample_n(1) # randomly sample point from each vell
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned) # run model
  predicted[[i]] <- avg_predictions(lm.thinned, by=c("abslat", "urban2"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                       urban2=c("Natural", "Suburban", "Urban"))) # store predictions for plotting
  ggeffects.slopes[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"), test = NULL) # store slopes
  ggeffects.slopes.contrast[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2")) # store contrasts
}

predicted_df <- bind_rows(predicted)
write.csv(predicted_df, "supplement_figs/90.coverage.thinned.results.csv")

ggeffects.slopes.df <- bind_rows(ggeffects.slopes)
write.csv(ggeffects.slopes.df, file="supplement_figs/90.coverage.thinned.slopes.csv")
ggslopes <- ggeffects.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), 
                                                                   conf.low=min(conf.low))
ggslopes


ggeffects.contrast_df <- bind_rows(ggeffects.slopes.contrast)
#write.csv(ggeffects.contrast_df, file="thinned_results/thinned_contrasts.csv")
contrast <- ggeffects.contrast_df %>% filter(urban2=="Suburban-Urban") %>% filter(p.value<0.05) # 1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Suburban") %>% filter(p.value<0.05) # 1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Urban") %>% filter(p.value<0.05) # 1000

######## Plot results
thinned.results.90 <- read.csv("supplement_figs/90.coverage.thinned.results.csv")
thinned.results.summary.90 <- thinned.results.90 %>% group_by(abslat, urban2) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))

thinned.plots.90 <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat_summary90, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary.90, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(thinned.results.summary.90, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,600))+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=18), 
        axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
# this is the plot with the 95% of the confidence intervals
thinned.plots.90
#ggsave(thinned.plots.90, file="supplement_figs/thinned.results.coverage90.png", height=6, width=7)




































##### 2) Different thresholds of modification
lowmod <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtLowThreshold.tif")
highmod <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtHighThreshold.tif")


## Do highmod first because I can just filter it from the modeling data (this is removing everything with more than values of 0.25 HI)
dat_highmod <- read.csv("modeling_data.csv")
# extract raster data
dat_highmod$urban <- as.data.frame(terra::extract(highmod, dat_highmod[,c(29:30)]))$SMOD_global
dat_highmod %>% group_by(urban) %>% count

# remove ones with NaN urbanization score
dat_highmod <- dat_highmod %>% filter(!is.na(urban)) # 42659 datapoints now

## Run model
high.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat_highmod)
anova(high.model)

square <- function(x){
  x^2
} 

## plot results
predicted.highmod<-avg_predictions(high.model, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))


GHMhighmod.plot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat_highmod, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.highmod, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.highmod, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
GHMhighmod.plot
## yeah pretty much the same
plot_slopes(high.model, variables="abslat", condition=c("urban2"))


### Thinned
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL)) # sample grid

# assign cell number to each point in my data
vect <- st_as_sf(dat_highmod, crs=st_crs(GHSL), coords=c("x","y"))
xy=st_coordinates(vect)
# get cell number that each point is in
dat_highmod$cell.subsample<-cellFromXY(sample.grid, xy)


square <- function(x){
  x^2
} # make function to square
predicted <- list()
ggeffects.slopes <- list()
ggeffects.slopes.contrast <- list()
set.seed(30)

for (i in 1:1000){
  dat.thinned <- dat_highmod %>% group_by(cell.subsample, urban2) %>% sample_n(1) # randomly sample point from each vell
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned) # run model
  predicted[[i]] <- avg_predictions(lm.thinned, by=c("abslat", "urban2"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                       urban2=c("Natural", "Suburban", "Urban"))) # store predictions for plotting
  ggeffects.slopes[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"), test = NULL) # store slopes
  ggeffects.slopes.contrast[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2")) # store contrasts
}

predicted_df <- bind_rows(predicted)
write.csv(predicted_df, "supplement_figs/high.mod.thinned.results.csv")

ggeffects.slopes.df <- bind_rows(ggeffects.slopes)
write.csv(ggeffects.slopes.df, file="supplement_figs/high.mod.thinned.slopes.csv")
ggslopes <- ggeffects.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), 
                                                                   conf.low=min(conf.low))
ggslopes


ggeffects.contrast_df <- bind_rows(ggeffects.slopes.contrast)
#write.csv(ggeffects.contrast_df, file="thinned_results/thinned_contrasts.csv")
contrast <- ggeffects.contrast_df %>% filter(urban2=="Urban-Suburban") %>% filter(p.value<0.05) # 1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Suburban") %>% filter(p.value<0.05) # 1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Urban-Natural") %>% filter(p.value<0.05) # 1000

######## Plot results
thinned.results.highmod <- read.csv("supplement_figs/high.mod.thinned.results.csv")
thinned.results.summary.highmod <- thinned.results.highmod %>% group_by(abslat, urban2) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))

thinned.plots.highmod <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat_highmod, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary.highmod, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(thinned.results.summary.highmod, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,600))+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=18), axis.title=element_blank())
# this is the plot with the 95% of the confidence intervals
thinned.plots.highmod
#ggsave(thinned.plots.highmod, file="supplement_figs/thinned.results.highmod.png", height=6, width=7)












### Now with lower threshold of human modification
summary.95 <- read.csv("global_richness_summary.csv") %>% filter(number_checklists >= 83) 

summary.95$x <- xFromCell(lowmod, summary.95$cell) # extract the coordinates from the cells
summary.95$y <- yFromCell(lowmod, summary.95$cell)

summary.95$urban <- as.data.frame(terra::extract(lowmod, summary.95[,c(21:22)]))$SMOD_global

summary.95 %>% group_by(urban) %>% summarise(n=n())

# remove ones with NaN urbanization score
summary.95_filt <- summary.95 %>% filter(!is.na(urban)) # 84,848 cells

## add all the other stuff in


# Turn into sf (spatial) object
summary.95_sf <- st_as_sf(summary.95_filt, coords=c("x", "y"), crs=st_crs(GHSL))

## Add elevation
# need to convert the points back to lat long
summary_latlong <- st_transform(summary.95_sf, crs=st_crs(4326))
latlong_df <- summary_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                         lat = sf::st_coordinates(.)[,2])

latlong_df <- get_elev_point(latlong_df[,c(23,24,2:21)], prj=crs(summary_latlong), src="aws", overwrite=TRUE) # extract elevations from amazon web services

dat_summary95 <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                                  y = sf::st_coordinates(.)[,2]))
#### assign hemisphere
dat_summary95$hemisphere <- "northern"
dat_summary95$hemisphere[dat_summary95$lat<0]<-"southern"

dat_summary95

## add precipitation
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
dat_summary95$precip <- as.data.frame(terra::extract(precip, dat_summary95[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12

dat_summary95$abslat <- abs(dat_summary95$lat)

## turn urban into 3 categories
dat_summary95 <- dat_summary95 %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat_summary95 %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat_summary95$urban2 <- as.factor(dat_summary95$urban2) # turn into factor for modelling

dat_summary95$urban2 <- factor(dat_summary95$urban2, levels = c("1", "2", "3"),
                               labels = c("Natural", "Suburban", "Urban")) # rename as natural, suburban, and urban

dat_lowmod <- dat_summary95

##### Model
low.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat_lowmod)

anova(low.model)

square <- function(x){
  x^2
} 

## plot results
predicted.lowmod<-avg_predictions(low.model, by=c("abslat", "urban2"), transform=square, 
                                   newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))


GHMlowmod.plot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat_lowmod, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.lowmod, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.lowmod, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
GHMlowmod.plot

plot_slopes(low.model, variables="abslat", condition=c("urban2"))



### Thinned
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL)) # sample grid

# assign cell number to each point in my data
vect <- st_as_sf(dat_lowmod, crs=st_crs(GHSL), coords=c("x","y"))
xy=st_coordinates(vect)
# get cell number that each point is in
dat_lowmod$cell.subsample<-cellFromXY(sample.grid, xy)


square <- function(x){
  x^2
} # make function to square
predicted <- list()
ggeffects.slopes <- list()
ggeffects.slopes.contrast <- list()
set.seed(30)

for (i in 1:1000){
  dat.thinned <- dat_lowmod %>% group_by(cell.subsample, urban2) %>% sample_n(1) # randomly sample point from each vell
  lm.thinned <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned) # run model
  predicted[[i]] <- avg_predictions(lm.thinned, by=c("abslat", "urban2"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                       urban2=c("Natural", "Suburban", "Urban"))) # store predictions for plotting
  ggeffects.slopes[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"), test = NULL) # store slopes
  ggeffects.slopes.contrast[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2")) # store contrasts
}

predicted_df <- bind_rows(predicted)
write.csv(predicted_df, "supplement_figs/low.mod.thinned.results.csv")

ggeffects.slopes.df <- bind_rows(ggeffects.slopes)
write.csv(ggeffects.slopes.df, file="supplement_figs/low.mod.thinned.slopes.csv")
ggslopes <- ggeffects.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), 
                                                                   conf.low=min(conf.low))
ggslopes


ggeffects.contrast_df <- bind_rows(ggeffects.slopes.contrast)
#write.csv(ggeffects.contrast_df, file="thinned_results/thinned_contrasts.csv")
contrast <- ggeffects.contrast_df %>% filter(urban2=="Suburban-Urban") %>% filter(p.value<0.05) # 999
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Suburban") %>% filter(p.value<0.05) # 1000
contrast <- ggeffects.contrast_df %>% filter(urban2=="Natural-Urban") %>% filter(p.value<0.05) # 1000

######## Plot results
thinned.results.lowmod <- read.csv("supplement_figs/low.mod.thinned.results.csv")
thinned.results.summary.lowmod <- thinned.results.lowmod %>% group_by(abslat, urban2) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))

thinned.plots.lowmod <- ggplot()+
  # geom_point(predicted.mean, mapping=aes(x=x, y=mean_x, color=group))+
  geom_point(dat_lowmod, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(thinned.results.summary.lowmod, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(thinned.results.summary.lowmod, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.5)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600), limits=c(0,600))+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=18), axis.title=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
# this is the plot with the 95% of the confidence intervals
thinned.plots.lowmod
#ggsave(thinned.plots.lowmod, file="supplement_figs/thinned.results.lowmod.png", height=6, width=7)


figure <- ggarrange(thinned.plots.98, thinned.plots.90, thinned.plots.highmod, thinned.plots.lowmod)
figure
figure2 <- annotate_figure(figure, left = textGrob("Species richness", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                           bottom = textGrob("Absolute latitude", gp = gpar(cex = 1.3)))

ggsave(figure2, file="supplement_figs/thinned.sensitivity.results.png", height=9, width=12)





