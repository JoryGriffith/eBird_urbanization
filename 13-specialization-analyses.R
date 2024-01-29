########### Script for analyzing whether specialist species are being disproportionately lost at lower latitudes
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)
library(marginaleffects)
library(ggpubr)
library(grid)
library(cowplot)
library(tidyterra)
library(ggeffects)
library(patchwork)
## First I want to thin the data so that for each cell, there is only one row for each species 

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
          "r4c2", "r4c4") # take out r4c1 and r4c4 because they are not in the final modelling data


model.data <- read.csv("modeling_data.csv") 
# also going to load this because I want to filter for cells that are in the final dataset to save space
unique(model.data$square)
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif")

datalist.years <- list()
datalist.names <- list()


for (i in 16:16){ # come back to 16, 18
  for (j in 5:6) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", names[i], "_", years[j], "_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% model.data$cell) # filter for cells that are in the final dataset
    dat.filt$SCIENTIFIC.NAME<- as.character(dat.filt$SCIENTIFIC.NAME)
    dat.filt$OBSERVATION.DATE<- as.character(dat.filt$OBSERVATION.DATE)
    dat.filt$OBSERVER.ID<- as.character(dat.filt$OBSERVER.ID)
    dat.filt$SAMPLING.EVENT.IDENTIFIER <- as.character(dat.filt$SAMPLING.EVENT.IDENTIFIER)
    dat.filt$OBSERVATION.COUNT <- as.character(dat.filt$OBSERVATION.COUNT)
    dat.filt$GROUP.IDENTIFIER <- as.character(dat.filt$GROUP.IDENTIFIER)
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  dat_uniquesp$x <- xFromCell(GHSL, dat_uniquesp$cell) # extract the coordinates from the cells
  dat_uniquesp$y <- yFromCell(GHSL, dat_uniquesp$cell)
  dat_uniquesp$square=names[i]
  datalist.names[[i]] <- dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(dat_uniquesp)
}

global_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(global_unique_sp$SCIENTIFIC.NAME)) # 10,723 species
#write.table(global_unique_sp, "global_unique_species.txt", row.names=FALSE)


###################################################################################################

GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif")
#global_uniquesp <- read.table("global_unique_sp.txt", header=TRUE) # this worked better than csv
length(unique(global_uniquesp$SCIENTIFIC.NAME))
#global_uniquesp <- global_uniquesp %>% na.omit() # remove random row with NA (not sure why that is there)

# add in lat long points to more easily bin by latitude
dat_latlong <- st_as_sf(global_uniquesp, coords=c("x", "y"), crs=st_crs(GHSL))
dat_latlong <- st_transform(dat_latlong, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- as.data.frame(dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                     lat = sf::st_coordinates(.)[,2]))

# bind this with data in other crs
global_uniquesp <- cbind(global_uniquesp, latlong_df[,5:6])

#### Extract urban scores
global_uniquesp$urban <- as.data.frame(terra::extract(GHSL, global_uniquesp[,c(3:4)]))$SMOD_global
test <- global_uniquesp %>% na.omit(urban) # there are no NAs in urban, this is good
# turn urban into 3 categories
global_uniquesp <- global_uniquesp %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), "natural", ifelse(urban==30, "urban", "suburban")))

global_uniquesp2 <- global_uniquesp %>% filter(lat <= 70 & lat >=-55) # filter for latitudes included in my analysis
length(unique(global_uniquesp$SCIENTIFIC.NAME)) 
length(unique(global_uniquesp2$SCIENTIFIC.NAME)) # only lost 3 species with the latitude cutoff

write.table(global_uniquesp2, "global_unique_species.txt", row.names=FALSE)


### Merge with trait data

############ Habitat data
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")

sp_habitat <- merge(global_uniquesp2, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
length(unique(sp_habitat$SCIENTIFIC.NAME)) # 8,496 species

# save data with habitat breadth
write.table(sp_habitat, "unique_sp_habitatbreadth.txt", row.names=F)

######## Diet data 
diet<- read.csv("/Volumes/Expansion/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
global_uniquesp2 <- read.table("global_unique_species.txt", header=T)
sp_diet <- merge(global_uniquesp2, diet[, c(9,21,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific") # merge with species data
length(unique(sp_diet$SCIENTIFIC.NAME)) # 6,902 species

# save data with diet
write.table(sp_diet, "unique_sp_dietspec.txt", row.names=F)



#############################################################














##### NEW WAY


########################################################


#### Trying new way of looking at specialization! (By cell)
global_uniquesp <- read.table("global_unique_species.txt", header=TRUE) %>% filter(lat <= 70 & lat >=-55)
# this is a list of species for each cell
# Merge with habitat data
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")
#habitat$habitat.scaled.before <- scales::rescale(habitat$Habitat_breadth_IUCN)

uniquesp_habitat <- merge(global_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")

habitat.species.summary <- uniquesp_habitat %>% group_by(cell, long, lat, urban2) %>% drop_na(Habitat_breadth_IUCN) %>% 
  summarise(mean.habitat=mean(Habitat_breadth_IUCN)) 
# scale habitat then subtract from 1
min.habitat <- min(habitat.species.summary$mean.habitat)
max.habitat <- max(habitat.species.summary$mean.habitat)

#habitat.species.summary$habitat.scaled <- scales::rescale(habitat.species.summary$mean.habitat)

habitat.species.summary$abslat <- abs(habitat.species.summary$lat)

# merge with modeling data to get precipitation and elevation
SR.dat <- read.csv("modeling_data.csv")
habitat.species.summary <- inner_join(habitat.species.summary, SR.dat[,c(1,2,3,21,26,34)], by="cell") %>% drop_na() 
habitat.species.summary$abslat <- abs(habitat.species.summary$lat)


# Run model
spec.mod1 <- lm(mean.habitat ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, data=habitat.species.summary)
summary(spec.mod1)
spec.mod2 <- lm(log(mean.habitat) ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, data=habitat.species.summary)
summary(spec.mod2)
anova(spec.mod2)

spec.mod3 <- lm(sqrt(mean.habitat) ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, data=habitat.species.summary)
summary(spec.mod3)

# better when log transformed
# adding in elevation, precipitation, hemisphere, and number of checklists all made the model better
AIC(spec.mod1, spec.mod2, spec.mod3)
# log transformed is the best
# the scaled habitat looks bad, something has gone wrong

#plot(spec.mod1) # doesn't look great
hist(spec.mod2$residuals)
# the model is right skewed

habitat.sf <- st_as_sf(habitat.species.summary, coords=c("long", "lat")) 
plot(gstat::variogram(residuals(spec.mod2, "pearson") ~
                                                1, data = habitat.sf, cutoff = 200))



exponent <- function(x){
  exp(x)
} 
summary(spec.mod2)
predicted.habitat <- avg_predictions(spec.mod2, by=c("abslat", "urban2"), transform=exponent,
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban")))

habitatLDG <- ggplot()+
  geom_point(habitat.species.summary, mapping=aes(x=abslat, y=mean.habitat, color=urban2), size=1, alpha=0.1)+
  geom_line(predicted.habitat, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.habitat, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_classic()+
  labs(y="Mean habitat breadth", x="Absolute latitude")+
  theme(legend.title=element_blank(), text=element_text(size=15), legend.position="none", axis.title.x=element_blank())
# will need to account for spatial autocorrelation
habitatLDG


avg_slopes(spec.mod2, variables="abslat", by="urban2")
hypothesis_test(spec.mod2, terms=c("abslat", "urban2"), scale="response")



### Plot latitudinal gradient with points colored by specialization
square <- function(x){
  x^2
} 
full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, habitat.species.summary)

predicted.full<-avg_predictions(full.model, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban")))


# want to make this over a more limited set of predictors and somehow make not wiggly
#marginal.full<-plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) # need to figure out how to back transform this
library(RColorBrewer)
mainLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(habitat.species.summary, mapping=aes(x=abslat, y=total_SR, shape=urban2, color=mean.habitat), size=1, alpha=0.8)+
  scale_shape_manual(values=c(3,17,19))+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, group=urban2, lty=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  labs(x="Absolute latitude", y="Species richness")+
#  scale_color_viridis_c(option = "magma")+ 
  scale_color_distiller(palette="Spectral")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mainLDGplot

library(spatialreg)
library(spdep)
# check for autocorrelation
habitat.samp <- habitat.species.summary[sample(nrow(habitat.species.summary), 5000), ]

spec.modsamp <- lm(sqrt(mean.habitat)~abslat * urban2 * hemisphere + elevation + precip + log(number_checklists), habitat.samp)
#library(sf)
habitat.sf <- st_as_sf(habitat.samp, coords=c("long", "lat")) 
habitat.nb <- dnearneigh(habitat.sf, d1=0, d2=200) # calculate distances
habitat.lw <- nb2listw(habitat.nb, style = "W", zero.policy = TRUE)
#library(mgcv)
spec.modsamp <- gls(log(mean.habitat)~abslat * urban2 * hemisphere + elevation + precip + log(number_checklists), habitat.samp)
plot(gstat::variogram(residuals(spec.modsamp, "pearson") ~
                        1, data = habitat.sf, cutoff = 200))
#
moran.samp <- lm.morantest(spec.modsamp, habitat.lw, zero.policy = T)
moran.samp

# significantly autocorrelated :(


#### Model habitat as a GAM
###################### Model habitat as a GAM
library(mgcv)
# model with regular latitude
habitat.species.summary$urban2 <- as.factor(habitat.species.summary$urban2)
mod.gam1 <- gam(mean.habitat ~ s(abslat, by=c(urban2)) + hemisphere + precip + log(number_checklists) + elevation, 
                data = habitat.species.summary, method="REML")
#mod.gam2 <- gam(mean.habitat ~ s(abslat, by=c(urban2), bs='ts') + hemisphere + precip + log(number_checklists) + elevation, data = habitat.species.summary)
AIC(mod.gam1, spec.mod2)
# the AIC for the GAM is lower
plot(ggeffects::ggpredict(mod.gam1, terms = c("abslat", "urban2")), facets = FALSE, add.data=TRUE, dot.alpha=0.2) 

#plot(ggeffects::ggpredict(mod.gam2, terms = c("abslat", "urban2")), facets = FALSE) 


habitat.sf <- st_as_sf(habitat.species.summary, coords=c("long", "lat")) 
# plot variogram
plot(gstat::variogram(residuals(spec.mod1, "response") ~
                        1, data = habitat.sf, cutoff = 200))
plot(gstat::variogram(residuals(beta.mod1.full, "response") ~
                        1, data = habitat.sf, cutoff = 200))
# the variogram of the linear model and the variogram of the beta mod look very similar

predicted.habitat.beta <- avg_predictions(beta.mod1.full, by=c("abslat", "urban2"),
                                          newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban")))










######################## Spatially thin models

#### Habitat
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL))

vect <- st_as_sf(habitat.species.summary, crs=st_crs(4326), coords=c("long","lat"))
vect2 <- st_transform(vect, crs=crs(GHSL))
xy=st_coordinates(vect2)
# get cell number that each point is in
habitat.species.summary$cell.subsample<-cellFromXY(sample.grid, xy)

# randomly sample one point within each cell
habitat.thinned <- habitat.species.summary %>% group_by(cell.subsample, urban2) %>% sample_n(1) 
test1 <- habitat.thinned %>% group_by(urban2) %>% summarise(mean = mean(mean.habitat))
filtered.out <- habitat.species.summary %>% filter(!cell %in% habitat.thinned$cell) # there are only 2,370 urban (so many)
test2 <- filtered.out %>% group_by(urban2) %>% summarise(mean = mean(mean.habitat)) # 11290 urban filtered out
# a lot of urban areas are filtered out and the mean habitat breadth is wider in the filtered out data



p1 <- ggplot(filtered.out, aes(x=abslat, y=mean.habitat, group=urban2))+
  geom_point(aes(color=urban2), size=1, alpha=0.3)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))

p2 <- ggplot(habitat.thinned, aes(x=abslat, y=mean.habitat, group=urban2))+
  geom_point(aes(color=urban2), size=1, alpha=0.3)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))

p1 | p2


habitat.thinned.mod1 <- lm(mean.habitat ~ abslat * urban2 * hemisphere + elevation + precip + log(number_checklists), data=habitat.thinned)
# whoa the second one is a lot better

habitat.thinned.sf <- st_as_sf(habitat.thinned, coords=c("long", "lat")) 
plot(gstat::variogram(residuals(habitat.thinned.mod1, "response") ~
                        1, data = habitat.thinned.sf, cutoff = 1000, width=20, alpha = c(0, 45, 90, 135)))


hist(residuals(habitat.thinned.mod1))
plot(habitat.thinned.mod1)
# plot thinned model
predicted.habitat.thin <- avg_predictions(habitat.thinned.mod1, by=c("abslat", "urban2"),
                                           newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban")))

habitatLDGthin <- ggplot()+
  geom_point(habitat.thinned, mapping=aes(x=abslat, y=mean.habitat, color=urban2), size=0.5, alpha=0.2)+
  geom_line(predicted.habitat.thin, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.habitat.thin, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_classic()+
  labs(y="Mean habitat specialization", x="Absolute latitude")+
  theme(legend.title=element_blank(), text=element_text(size=15), legend.position="none", axis.title.x=element_blank())
# will need to account for spatial autocorrelation
habitatLDGthin


## Look at GAM of thinned data
library(mgcv)
# model with regular latitude
habitat.thinned$urban2 <- as.factor(habitat.thinned$urban2)
mod.gam1 <- gam(mean.habitat ~ s(abslat, by=c(urban2)) + hemisphere + precip + log(number_checklists) + elevation, data = habitat.thinned, method="REML")

plot(ggeffects::ggpredict(mod.gam1, terms = c("abslat", "urban2")), facets = FALSE, add.data=TRUE, alpha=0.2) 
AIC(mod.gam1, habitat.thinned.mod1)
# linear model is better

# try with mean habitat
habitat.thinned$urban2 <- as.factor(habitat.thinned$urban2)
mod.gam1 <- gam(mean.habitat ~ s(abslat, by=c(urban2)) + hemisphere + precip + log(number_checklists) + elevation, data = habitat.thinned)

# try beta regression


plot(ggeffects::ggpredict(habitat.thinned.mod1, terms = c("abslat", "urban2")), facets = FALSE) 
# very similar output between the linear model and the beta regression model
summary(beta.mod1) # the R2 is much better than the linear model
summary(habitat.thinned.mod1)
AIC(beta.mod3, beta.mod4) # mod 3 is the best


predicted.full<-avg_predictions(beta.mod3, by=c("abslat", "urban2"),
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban")))


ggplot()+
#  geom_point(total.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  geom_line(predicted.habitat.beta, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5, lty=2)+
  geom_ribbon(predicted.habitat.beta, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  theme_classic()
# the dotted line is the line with the full data and the other line is the thinned data





AIC(habitat.thinned.mod1, beta.mod1) # the beta regression is better (even then the log transformed)




####### Start iteratively thinning
exponent <- function(x){
  exp(x)
}
predicted <- list()
means <- list()
ggeffects.slopes <- list()
ggeffects.slopes.contrast <- list()
ggeffects.slopes.contrast.hemisphere <- list()
set.seed(100)

for (i in 1:1000){
  dat.thinned <- habitat.species.summary %>% group_by(cell.subsample, urban2) %>% sample_n(1) 
  lm.thinned <- lm(log(mean.habitat) ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned)
  # dat.thinned.sf <- st_as_sf(dat.thinned, coords=c("long", "lat")) 
  #dat.thinned.nb <- dnearneigh(dat.thinned.sf, d1=0, d2=200) # calculate distances
  #  dat.thinned.lw <- nb2listw(dat.thinned.nb, style = "W", zero.policy = TRUE) # turn into weighted list
  # moran <- lm.morantest(lm.thinned, dat.thinned.lw, zero.policy = T)
  predicted[[i]] <- avg_predictions(lm.thinned, by=c("abslat", "urban2"), transform=exponent, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban"))) # store predictions for plotting
  means[[i]] <- marginal_means(lm.thinned, variables=c("abslat", "urban2"), transform=exponent)
  #  means[[i]] <- emmeans(lm.thinned, specs="urban2")
  ggeffects.slopes[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"), test = NULL) # see differences in slopes (using ggeffects)
  ggeffects.slopes.contrast[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2"))
  ggeffects.slopes.contrast.hemisphere[[i]] <- hypothesis_test(lm.thinned, c("abslat", "urban2", "hemisphere"))
}

predicted_df <- bind_rows(predicted)
write.csv(predicted_df, "thinned_results/thinned.habitat.specialization.csv")

ggeffects.slopes.df <- bind_rows(ggeffects.slopes)
ggslopes <- ggeffects.slopes.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), conf.low=min(conf.low))
ggslopes
write.csv(ggeffects.slopes.df, "thinned_results/thinned.habitat.specialization.slopes.csv")
# urban steeper than natural, natural second steepest, suburban the least steep


ggeffects.slopes.contrast.df <- bind_rows(ggeffects.slopes.contrast)
contrast <- ggeffects.slopes.contrast.df %>% filter(urban2=="urban-suburban") %>% filter(p.value<0.05) # 1000
contrast <- ggeffects.slopes.contrast.df %>% filter(urban2=="natural-suburban") %>% filter(p.value<0.05) #1000
contrast <- ggeffects.slopes.contrast.df %>% filter(urban2=="natural-urban") %>% filter(p.value<0.05) # 74
# all significantly different slopes

ggeffects.contrast.hemisphere_df <- bind_rows(ggeffects.slopes.contrast.hemisphere)

contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="natural-natural", hemisphere=="northern-southern") %>% filter(p.value<0.05) # 1000
# steeper in N hemsiphere
contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="urban-urban", hemisphere=="northern-southern") %>% filter(p.value<0.05) # 991
contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="suburban-suburban", hemisphere=="northern-southern") %>% filter(p.value<0.05) # 0

contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="urban-suburban", hemisphere=="southern-southern") %>% filter(p.value<0.05) # 54
# 718
contrast <- ggeffects.contrast.hemisphere_df %>% filter(urban2=="urban-suburban", hemisphere=="northern-northern") %>% filter(p.value<0.05) # 1000



########## Plot
predicted_df <- read.csv("thinned_results/thinned.habitat.specialization.csv")
habitat.thinned.results.summary <- predicted_df %>% group_by(abslat, urban2) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))


thinned.habitatLDG <- ggplot()+
  geom_point(habitat.species.summary, mapping=aes(x=abslat, y=mean.habitat, color=urban2), size=0.25, alpha=0.1)+
  geom_line(habitat.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(habitat.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_classic()+
  labs(y="Mean habitat breadth", x="Absolute latitude")+
  theme(legend.title=element_blank(), text=element_text(size=15), legend.position="none", axis.title.x=element_blank())
# will need to account for spatial autocorrelation
thinned.habitatLDG


11.74-9.56 # 2.18 (diff between habitat breadth at equator)

18.82-15.78 # 3.04 (diff between habitat breadth at equator)




###### Look at relationship between number of cells in a larger grid cell and specialization (are larger cities more specialized?)

cells.no <- habitat.species.summary %>% group_by(cell.subsample, urban2) %>% summarise(n=n(), mean.habitat=mean(mean.habitat)) %>% filter(urban2=="urban")

ggplot(cells.no, aes(x=n, y=mean.habitat))+
  geom_point()+
  geom_smooth(method="lm")



















########## Diet

##### Now try with diet
global_uniquesp <- read.table("global_unique_species.txt", header=TRUE) %>% filter(lat <= 70 & lat >=-55)
diet<- read.csv("/Volumes/Expansion/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
uniquesp_diet <- merge(global_uniquesp, diet[, c(9,21,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific")
range(uniquesp_diet$gini.index) # 0.35 to 1, pretty goof range
uniquesp_diet$gini.flipped <- 1-uniquesp_diet$gini.index
range(uniquesp_diet$gini.flipped) #0 to 0.644
diet.species.summary <- uniquesp_diet %>% group_by(cell, long, lat, urban2) %>% drop_na(gini.flipped) %>% 
  summarise(mean.diet.flipped=mean(gini.flipped), mean.diet=mean(gini.index))
range(diet.species.summary$mean.diet.flipped) # 0.017 to 0.245
range(diet.species.summary$mean.diet)  # 0.75 to 0.98 
# this is the same difference
diet.species.summary$abslat <- abs(diet.species.summary$lat)
SR.dat <- read.csv("modeling_data.csv")
diet.species.summary <- inner_join(diet.species.summary, SR.dat[,c(1,2,3,21,26,34)], by="cell") %>% drop_na()


diet.mod1 <- lm(mean.diet.flipped ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, data=diet.species.summary)
# R2 is 0.39, pretty good!!
# AIC is the same when it is flipped
# model without log transformation is better
#AIC(diet.mod1, diet.mod2)

#plot(diet.mod1)
summary(diet.mod1)
hist(diet.mod1$residuals) # they honestly look pretty good


predicted.diet <- avg_predictions(diet.mod1, by=c("abslat", "urban2"),
                                     newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                        urban2=c("natural", "suburban", "urban")))

dietLDG <- ggplot()+
  geom_point(diet.species.summary, mapping=aes(x=abslat, y=mean.diet.flipped, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.diet, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.diet, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_classic()+
  labs(y="Mean diet specialization", x="Absolute latitude")+
  theme(legend.title=element_blank(), text=element_text(size=15), legend.position="none", axis.title.x=element_blank())
dietLDG
#specializationLDG <- habitatLDG / dietLDG
# models are the same just flipped!
#ggsave(specializationLDG, file="specializationLDG.png", height=7, width=5)

hypothesis_test(diet.mod1, terms=c("abslat", "urban2"), scale="response")
avg_slopes(diet.mod1, variables="abslat", by="urban2")


#### Plot full LDG colored by diet specialization

full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, diet.species.summary)

predicted.full<-avg_predictions(full.model, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban")))


# want to make this over a more limited set of predictors and somehow make not wiggly
#marginal.full<-plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) # need to figure out how to back transform this
library(RColorBrewer)
mainLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(diet.species.summary, mapping=aes(x=abslat, y=total_SR, shape=urban2, color=mean.diet.flipped), size=1, alpha=0.8)+
  scale_shape_manual(values=c(3,17,19))+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, group=urban2, lty=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  labs(x="Absolute latitude", y="Species richness")+
  #  scale_color_viridis_c(option = "magma")+ 
  scale_color_distiller(palette="Spectral")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mainLDGplot




### Run GAM with diet specialization
library(mgcv)
diet.species.summary$urban2 <- as.factor(diet.species.summary$urban2)
mod.gam1 <- gam(mean.diet.flipped ~ s(abslat, by=c(urban2)) + hemisphere + precip + log(number_checklists) + elevation, data = diet.species.summary)

plot(ggeffects::ggpredict(mod.gam1, terms = c("abslat", "urban2")), facets = FALSE, add.data=TRUE, alpha=0.2)
AIC(diet.mod1, mod.gam1)
# linear model is better

library(RColorBrewer)
test <- diet.species.summary %>% group_by(cell.subsample) %>% mutate(n=n())
ggplot(test, aes(x=abslat, y=mean.diet, color=n, shape=urban2))+
  geom_point(alpha=0.3)+
  scale_color_distiller(palette="Spectral")

### Try beta regression
library(betareg)
beta.mod1 <- betareg(mean.diet.flipped ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, data=diet.species.summary, link="log")
#summary(beta.mod1) # everything is significant
plot(beta.mod1, which = 5, type = "pearson")
AIC(diet.mod1, mod.gam1, beta.mod1)
# all of the betareg models are worse than the regular linear model




#####################################################





#############################
## Thinning
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_global.tif")
spat.extent <- ext(GHSL)
sample.grid <- rast(resolution=c(10000, 10000), extent = spat.extent, crs=crs(GHSL))

vect <- st_as_sf(diet.species.summary, crs=st_crs(4326), coords=c("long","lat"))
vect2 <- st_transform(vect, crs=crs(GHSL))
xy=st_coordinates(vect2)
# get cell number that each point is in
diet.species.summary$cell.subsample<-cellFromXY(sample.grid, xy)

# randomly sample one point within each cell
diet.thinned <- diet.species.summary %>% group_by(cell.subsample, urban2) %>% sample_n(1) 
test1 <- diet.thinned %>% group_by(urban2) %>% summarise(mean = mean(mean.diet.flipped))
filtered.out <- diet.species.summary %>% filter(!cell %in% diet.thinned$cell) 
test2 <- filtered.out %>% group_by(urban2) %>% summarise(mean = mean(mean.diet.flipped)) 

p1 <- ggplot(filtered.out, aes(x=abslat, y=mean.diet, group=urban2))+
  geom_point(aes(color=urban2), size=1, alpha=0.3)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))

p2 <- ggplot(diet.thinned, aes(x=abslat, y=mean.diet, group=urban2))+
  geom_point(aes(color=urban2), size=1, alpha=0.3)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))

p1 | p2
# ok they are less different now!


lm.thinned.diet <- lm(mean.diet.flipped ~ abslat * urban2 * hemisphere +
                        precip + log(number_checklists) + elevation, diet.thinned)


hist(residuals(lm.thinned.diet))


## GAM of thinned
diet.thinned$urban2 <- as.factor(diet.thinned$urban2)
mod.gam1 <- gam(mean.diet.flipped ~ s(abslat, by=c(urban2)) + hemisphere + precip + log(number_checklists) + elevation, data = diet.thinned)
summary(mod.gam1)
plot(ggeffects::ggpredict(mod.gam1, terms = c("abslat", "urban2")), facets = FALSE, add.data=TRUE, alpha=0.2) 


##### Beta regression on thinned
library(betareg)
beta.mod1 <- betareg(mean.diet.flipped ~ abslat * urban2 * hemisphere + precip + log(number_checklists) + elevation, data=diet.thinned, link="log")

summary(beta.mod1) # everything is significant
plot(beta.mod1, which = 5, type = "pearson")
hist(residuals(beta.mod1))
# log link is the best
plot(ggeffects::ggpredict(beta.mod1, terms = c("abslat", "urban2")), facets = FALSE, add.data=TRUE, alpha=0.2) 
plot(ggeffects::ggpredict(lm.thinned.diet, terms = c("abslat", "urban2")), facets = FALSE, add.data=TRUE, alpha=0.2)
AIC(lm.thinned.diet, beta.mod1)
# linear model is better


# plot randomly sampled
predicted.diet.thin <- avg_predictions(beta.mod1, by=c("abslat", "urban2"),
                                  newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), 
                                                     urban2=c("natural", "suburban", "urban")))

dietLDGthinned <- ggplot()+
  geom_point(diet.species.summary, mapping=aes(x=abslat, y=mean.diet.flipped, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.diet.thin, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.diet.thin, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  geom_line(predicted.diet, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5, lty=2)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_classic()+
  labs(y="Mean diet breadth", x="Absolute latitude")+
  theme(legend.title=element_blank(), text=element_text(size=15), legend.position="none", axis.title.x=element_blank())

dietLDGthinned
plot(lm.thinned.diet)


####### Start thinning
predicted.diet <- list()
ggeffects.slopes.diet <- list()
ggeffects.slopes.contrast.diet <- list()
ggeffects.slopes.contrast.hemisphere.diet <- list()
set.seed(100)

for (i in 1:1000){
  dat.thinned.diet <- diet.species.summary %>% group_by(cell.subsample, urban2) %>% sample_n(1) 
  lm.thinned.diet <- lm(mean.diet.flipped ~ abslat * urban2 * hemisphere +
                     precip + log(number_checklists) + elevation, dat.thinned.diet)
  # dat.thinned.sf <- st_as_sf(dat.thinned, coords=c("long", "lat")) 
  #dat.thinned.nb <- dnearneigh(dat.thinned.sf, d1=0, d2=200) # calculate distances
  #  dat.thinned.lw <- nb2listw(dat.thinned.nb, style = "W", zero.policy = TRUE) # turn into weighted list
  # moran <- lm.morantest(lm.thinned, dat.thinned.lw, zero.policy = T)
  predicted.diet[[i]] <- avg_predictions(lm.thinned.diet, by=c("abslat", "urban2"),  
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("natural", "suburban", "urban"))) # store predictions for plotting
  ggeffects.slopes.diet[[i]] <- hypothesis_test(lm.thinned.diet, c("abslat", "urban2"), test = NULL) # see differences in slopes (using ggeffects)
  ggeffects.slopes.contrast.diet[[i]] <- hypothesis_test(lm.thinned.diet, c("abslat", "urban2"))
  ggeffects.slopes.contrast.hemisphere.diet[[i]] <- hypothesis_test(lm.thinned.diet, c("abslat", "urban2", "hemisphere"))
}

predicted_diet_df <- bind_rows(predicted.diet)
write.csv(predicted_diet_df, "thinned_results/thinned.diet.specialization.csv")

ggeffects.slopes.diet.df <- bind_rows(ggeffects.slopes.diet)
ggslopes.diet <- ggeffects.slopes.diet.df %>% group_by(urban2) %>% summarise(mean=mean(Slope), conf.high = max(conf.high), conf.low=min(conf.low))
ggslopes.diet
write.csv(ggeffects.slopes.diet.df, "thinned_results/thinned.diet.specialization.slopes.csv")
# urban steeper than natural, natural second steepest, suburban the least steep


ggeffects.slopes.contrast.diet.df <- bind_rows(ggeffects.slopes.contrast.diet)
contrast <- ggeffects.slopes.contrast.diet.df %>% filter(urban2=="urban-suburban") %>% filter(p.value<0.05) # 0
contrast <- ggeffects.slopes.contrast.diet.df %>% filter(urban2=="natural-suburban") %>% filter(p.value<0.05) #1000
contrast <- ggeffects.slopes.contrast.diet.df %>% filter(urban2=="natural-urban") %>% filter(p.value<0.05) # 1000
# all significantly different slopes

ggeffects.contrast.hemisphere.diet_df <- bind_rows(ggeffects.slopes.contrast.hemisphere.diet)

contrast <- ggeffects.contrast.hemisphere.diet_df %>% filter(urban2=="natural-natural", hemisphere=="northern-southern") %>% filter(p.value<0.05) # 1000
# steeper in N hemsiphere
contrast <- ggeffects.contrast.hemisphere.diet_df %>% filter(urban2=="urban-urban", hemisphere=="northern-southern") %>% filter(p.value<0.05) # 2
contrast <- ggeffects.contrast.hemisphere.diet_df %>% filter(urban2=="suburban-suburban", hemisphere=="northern-southern") %>% filter(p.value<0.05) # 757

contrast <- ggeffects.contrast.hemisphere.diet_df %>% filter(urban2=="urban-suburban", hemisphere=="southern-southern") %>% filter(p.value<0.05) # 0
# 718
contrast <- ggeffects.contrast.hemisphere.diet_df %>% filter(urban2=="urban-suburban", hemisphere=="northern-northern") %>% filter(p.value<0.05) # 0


########## Plot
predicted_diet_df <- read.csv("thinned_results/thinned.diet.specialization.csv")
diet.thinned.results.summary <- predicted_diet_df %>% group_by(abslat, urban2) %>% 
  summarise(mean_x=mean(estimate), max.conf.high = max(conf.high), min.conf.low = min(conf.low))


thinned.dietLDG <- ggplot()+
  geom_point(diet.species.summary, mapping=aes(x=abslat, y=mean.diet.flipped, color=urban2), size=0.25, alpha=0.1)+
  geom_line(diet.thinned.results.summary, mapping=aes(x=abslat, y=mean_x, color=urban2), lwd=1.5)+
  geom_ribbon(diet.thinned.results.summary, mapping=aes(x=abslat, ymax=max.conf.high, ymin=min.conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  theme_classic()+
  labs(y="Mean diet breadth", x="Absolute latitude")+
  theme(legend.title=element_blank(), text=element_text(size=15), legend.position="none", axis.title.x=element_blank())
# will need to account for spatial autocorrelation
thinned.dietLDG


thinned.specialization.plots <- thinned.habitatLDG / thinned.dietLDG
thinned.specialization.plots
ggsave(thinned.specialization.plots, file="thinned.specializationLDG.png", height=7, width=5)


specializationLDG | thinned.specialization.plots



###### Look at relationship between number of cells in a larger grid cell and specialization (are larger cities more specialized?)

cells.no <- diet.species.summary %>% group_by(cell.subsample, urban2) %>% summarise(n=n(), mean.diet=mean(mean.diet)) %>% filter(urban2=="urban")

ggplot(cells.no, aes(x=log(n), y=log(mean.diet)))+
  geom_point()+
  geom_smooth(method="lm")














######## Make density distributions of species at different binned latitudes ########

## habitat
global_uniquesp <- read.table("global_unique_species.txt", header=TRUE) %>% filter(lat <= 70 & lat >=-55)
# this is a list of species for each cell
# Merge with habitat data
habitat <- read.csv("/Volumes/Backup/eBird/Traits/habitat_breadth.csv")
diet<- read.csv("/Volumes/Backup/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
#habitat$habitat.scaled.before <- scales::rescale(habitat$Habitat_breadth_IUCN)

uniquesp_habitat <- merge(global_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
uniquesp_habitat$abslat <- abs(uniquesp_habitat$lat)

# bin by latitude zones
uniquesp_habitat <- uniquesp_habitat %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 70), 
                                                                labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))
# make list of unique species in each zone bin
unique_zone_bin <- uniquesp_habitat %>% distinct(zone_bin, urban2, SCIENTIFIC.NAME)
# merge with habitat data
unique_zone_bin <- merge(unique_zone_bin, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial", all.x=TRUE)
?merge
# merge with diet data
unique_zone_bin <- merge(unique_zone_bin, diet[, c(9,21,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific", all.x=TRUE)
unique_zone_bin$gini.flipped <- 1-unique_zone_bin$gini.index

### Plot density plots
## Habitat
ggplot(unique_zone_bin, aes(x = Habitat_breadth_IUCN, y=after_stat(count), fill=urban2)) +
  geom_density(alpha=0.6) +
  facet_wrap(~zone_bin)
# this makes much more sense

## Diet
ggplot(unique_zone_bin, aes(x = gini.flipped, y=after_stat(count), fill=urban2)) +
  geom_density(alpha=0.6) +
  facet_wrap(~zone_bin)




#### also try binning by smaller zones
uniquesp_habitat <- uniquesp_habitat %>% mutate(small_bin = cut(abslat, breaks=c(0, 10, 20, 30, 40, 50, 60, 70)))
unique_zone_small <- uniquesp_habitat %>% distinct(small_bin, urban2, SCIENTIFIC.NAME)
# merge with habitat data
unique_zone_small <- merge(unique_zone_small, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial", all.x=TRUE)

# merge with diet data
unique_zone_small <- merge(unique_zone_small, diet[, c(9,21,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific", all.x=TRUE)
unique_zone_small$gini.flipped <- 1-unique_zone_small$gini.index


### Plot density plots
## Habitat
ggplot(unique_zone_small, aes(x = log(Habitat_breadth_IUCN), y=after_stat(count), fill=urban2)) +
  geom_density(alpha=0.6) +
  facet_wrap(~small_bin)
# this makes much more sense

## Diet
ggplot(unique_zone_small, aes(x = gini.flipped, y=after_stat(count), fill=urban2)) +
  geom_density(alpha=0.6) +
  facet_wrap(~small_bin)









########################################################
########### OLD WAY #######################
#######################################






## Habitat data, filter out NAs for habitat
sp_habitat <- read.table("unique_sp_habitatbreadth.txt", header=T) %>% filter(!is.na(Habitat_breadth_IUCN)) # some species have NA values for habitat breadth
length(unique(sp_habitat$SCIENTIFIC.NAME)) # 8367

# see if there are more specialists at low latitudes
# boxplot of habitat specialization
hist(log(sp_habitat$Habitat_breadth_IUCN)) # definitely looks pretty log normal
hist(sp_habitat$Habitat_breadth_IUCN)
sp_habitat$abslat <- abs(sp_habitat$lat) # add in absolute latitude


######### Bin by biogeographical zone
sp_habitat <- sp_habitat %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 70), labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))

# run anova on the raw habitat breadth with these larger zones
habitat.aov3 <- aov(Habitat_breadth_IUCN ~ zone_bin * urban2, data = sp_habitat)
summary(habitat.aov3) # significant interaction

# look at just difference in specialization with latitude
marginal_means(habitat.aov3, variables=c("zone_bin"))

means.habitat3 <- marginal_means(habitat.aov3, variables=c("zone_bin", "urban2"))
# the difference is way larger in the tropics!



birds_zones <- sp_habitat %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% 
  pivot_wider(names_from="urban2", values_from="n")  

birds_zones <- birds_zones %>% replace(is.na(.), 0)

birds_zones$category <- NA

for (i in 1:nrow(birds_zones)){
  if (birds_zones$natural[i] > 0 & birds_zones$urban[i] > 0) {
    birds_zones$category[i] <- "in.urban"
  }
    else if (birds_zones$natural[i] == 0 & birds_zones$urban[i] > 0){
     birds_zones$category[i] <- "urban.only"
    }
  else if (birds_zones$natural[i] > 0 & birds_zones$urban[i] == 0) {
    birds_zones$category[i] <- "natural.only"
  }
  
}

# look at birds in different categories and see where they are found
urban.only <- birds_zones %>% filter(category=="urban.only") # remove urban only birds?
urban.only2 <- urban.only %>% pivot_wider(names_from="zone_bin", values_from="urban")


### denisty plot of trait values
ggplot(birds_zones, aes(x = log(Habitat_breadth_IUCN), y=after_stat(count), fill=category)) +
  geom_density(alpha=0.6) +
  facet_wrap(~zone_bin)





#both <- birds_zones %>% filter(category=="both")
#plot(both$natural~both$urban) # definitely a positive correlation
#cor(both$natural,both$urban) # 0.935

#length(unique(urb.only$SCIENTIFIC.NAME)) # 319 species that are urban only

# with suburban included
#for (i in 1:nrow(birds_zones)){
#  if (birds_zones$natural[i] >= 0 & birds_zones$suburban[i] >= 0 & birds_zones$urban[i] > 0) {
#    birds_zones$category[i] <- "in.urban"
#  }
#  else if (birds_zones$natural[i] >= 0 & birds_zones$suburban[i] > 0 & birds_zones$urban[i] == 0){
#    birds_zones$category[i] <- "in.suburban"
#  }
#  else if (birds_zones$natural[i] > 0 & birds_zones$suburban[i] == 0 & birds_zones$urban[i] == 0) {
#    birds_zones$category[i] <- "natural.only"
#  }
#  
#}



habitat.aov4 <- aov(log(Habitat_breadth_IUCN) ~ zone_bin * category, data = birds_zones)
summary(habitat.aov4)
# need to change transformation to exponent!
exponent <- function(x){
  exp(x)
} 
means.habitat4 <- marginal_means(habitat.aov4, variables=c("zone_bin", "category"), cross=TRUE, transform=exponent)
?marginal_means
#emmeans(habitat.aov4, pairwise~zone_bin, by="category") # all means are significantly different from one another across zones

means.habitat4
# Plot of richness of different categories

# this is exactly what I was trying to show! Didn't know it would be so clear
# when I add in urban only it doesn't really change anything because there are so few species that are urban only

# Make ridgeline plot
# Jenn asked me to do this but it is not telling that much to be honest
library(ggridges)
ggplot(birds_zones, aes(y = zone_bin, x = Habitat_breadth_IUCN, fill = category)) +
  geom_density_ridges() +
  theme_ridges() 

ggplot(birds_zones, aes(y = zone_bin, x = log(Habitat_breadth_IUCN), fill = category)) +
  geom_violin() +
  theme_ridges() 

# Stacked bar plot of birds in urban and not in urban
birds_zones %>% group_by(zone_bin) %>% count()
richness_category <- birds_zones %>% group_by(zone_bin, category) %>% count()

habitat_bar <-
  richness_category %>% 
  mutate(category = factor(category, levels = c('natural.only', 'in.urban', 'urban.only'), ordered = TRUE)) %>%
  ggplot(aes(fill=category, y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('In natural', 'In Urban', 'In urban only'), values=c("deepskyblue3", "grey30"))+
  labs(y="Number of Species")+
  coord_flip()+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title=element_blank(),legend.position = c(.95, .75),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
habitat_bar
#  theme(axis.title.x=element_blank(), legend.position="none")

### Plot emmeans using ggplot

habitat_point <- 
  ggplot(means.habitat4, aes(y=zone_bin, x=estimate, group=category, color=category))+
  geom_point(size=2)+
  geom_line(size=0.5)+
  scale_color_manual(values=c("grey30", "deepskyblue3"))+
  #  scale_y_reverse()+
  labs(x="Habitat breadth")+
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=0.1)+
  # annotate("text", x=0.7, y=2.5, label="Generalist", angle=90)+
  #  annotate("text", x=0.7, y=1.4, label="Specialist", angle=90)+
  # annotate("segment", x = 0.7, y = 2.75, xend = 0.7, yend = 2.9, size=0.6,
  #         arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
  #annotate("segment", x = 0.7, y = 1, xend = 0.7, yend = 1.15, size=0.5,
  #        arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="first"))+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(axis.title.y=element_blank(), legend.position="none") 
habitat_point
# theme(axis.title.x=element_blank(), legend.title=element_blank(),legend.position = c(.95, .1),
#      legend.justification = c("right", "bottom"),
#     legend.box.just = "right",
#    legend.margin = margin(6, 6, 6, 6))

habitat_plot <- ggarrange(habitat_bar, habitat_point, ncol=1)
habitat_plot


#ggsave(habitat_plot, file="pecialistHabitatResults.png", height=4, width=8)
# Try it as an inset




## Try making a plot using EulerR
birds_zones %>% group_by(zone_bin) %>% count()
richness_category <- birds_zones %>% group_by(zone_bin, category) %>% count()
# put it into a dataframe that eulerr will understand
library(eulerr)
euler_df <- c(Ntrop=2934, Utrop=113, "Ntrop&Utrop"=4251, Nsub=939, Usub=77, "Nsub&Usub"=2918, 
              Ntemp=383, Utemp=111, "Ntemp&Utemp"=1537, Npol=217, Upol=36, "Npol&Upol"=690) # tropical

euler_trop <- euler(euler_df, shape="circle")
plot(euler_trop)







## Plt distribution of some of these individually
# seen in the most natural cells but no urban or suburban - Cyrtonyx montezumae
global_uniquesp <- read.table("global_unique_species.txt", header=TRUE)
C.montezumae <- global_uniquesp %>% filter(SCIENTIFIC.NAME=="Cyrtonyx montezumae")

library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf", country=c("United States of America", "mexico"))
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filt_3cat.tif")

ggplot(data=world)+
  geom_sf() +
  geom_spatraster(data=GHSL$SMOD_global)+
  geom_point(data=C.montezumae, aes(x=long, y=lat), size=0.1, color="white") +
  coord_sf(expand = FALSE, xlim=c(-120, -80), ylim=c(10,40)) +
  labs(x="Longitude", y="Latitude")+
  theme_bw()

###########################################


### Diet specialization
sp_diet <- read.table("unique_sp_dietspec.txt", header=T) %>% filter(!is.na(gini.index))

length(unique(sp_diet$SCIENTIFIC.NAME)) # 
# see if there are more specialists at low latitudes


sp_diet$abslat <- abs(sp_diet$lat)
# boxplot of specialization


####### Bin by geographical zone

sp_diet <- sp_diet %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 70), labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))


# run anova on the raw habitat breadth with these larger zones
diet.aov3 <- aov(gini.index ~ zone_bin * urban2, data = sp_diet)
summary(diet.aov3) # significant interaction

means.diet3 <- marginal_means(diet.aov3, variables=c("zone_bin", "urban2"), cross=TRUE)
# the difference is way larger in the tropics!

diet_zones <- sp_diet %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, gini.index, Diet.5Cat) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>%
  pivot_wider(names_from="urban2", values_from="n")  

diet_zones <- diet_zones %>% replace(is.na(.), 0)

diet_zones$category <- NA

for (i in 1:nrow(diet_zones)){
  if (diet_zones$natural[i] >= 0  & diet_zones$urban[i] > 0) {
    diet_zones$category[i] <- "both"
  }
  # else if (diet_zones$natural[i] == 0 & diet_zones$urban[i] > 0) {
  # diet_zones$category[i] <- "urban.only"
  #}
  else if (diet_zones$natural[i] > 0 & diet_zones$urban[i] == 0) {
    diet_zones$category[i] <- "natural.only"
  }
}




ggplot(diet_zones, aes(x = log(gini.index), y=after_stat(count), fill=category)) +
  geom_density(alpha=0.6) +
  facet_wrap(~zone_bin)








# including suburban
#for (i in 1:nrow(diet_zones)){
#  if (diet_zones$natural[i] >= 0 & diet_zones$suburban[i] >= 0 & diet_zones$urban[i] > 0) {
#   diet_zones$category[i] <- "in.urban"
# }
#  else if (diet_zones$natural[i] >= 0 & diet_zones$suburban[i] > 0 & diet_zones$urban[i] == 0) {
#   diet_zones$category[i] <- "in.suburban"
# }

#else if (diet_zones$natural[i] > 0 & diet_zones$suburban[i] == 0 & diet_zones$urban[i] == 0) {
#   diet_zones$category[i] <- "natural.only"
#  }
#}


unique(diet_zones$zone_bin)
diet.aov4 <- aov(gini.index ~ zone_bin * category, data = diet_zones)
summary(diet.aov4)

means.diet4 <- marginal_means(diet.aov4, variables=c("zone_bin", "category"), cross=TRUE)
#emmeans(diet.aov4, pairwise~zone_bin, by="category")

##### Stacked bar plot of overall species loss
richness_category <- diet_zones %>% group_by(zone_bin, category) %>% count()
ggplot(richness_category, aes(fill=category, y=n, x=zone_bin)) + 
  geom_bar(position="stack", stat="identity") + # higher proportion of species being lost at the equator
  labs(y="richness")

library(ggridges)
ggplot(diet_zones, aes(y = zone_bin, x = gini.index, fill = category)) +
  geom_density_ridges() +
  theme_ridges() 

ggplot(diet_zones, aes(y = zone_bin, x = gini.index, fill = category)) +
  geom_violin() +
  theme_ridges()

##########



richness_category <- diet_zones %>% group_by(zone_bin, category) %>% count()

diet_bar <- richness_category %>% 
  mutate(category = factor(category, levels = c('natural.only', 'both', 'urban.only'), ordered = TRUE)) %>%
  ggplot(aes(fill=category, y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('In natural', 'In urban', 'In urban only'), values=c("deepskyblue3", "grey30"))+
  geom_bar(position="stack", stat="identity")+
  labs(y="Number of Species")+
  #  coord_flip()+
  theme_classic()+
  #  theme(axis.title.x=element_blank(), legend.position="none")
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title=element_blank(),legend.position = c(.95, .75),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
diet_bar



diet_point <- means.diet4 %>% 
  ggplot(aes(y=zone_bin, x=estimate, group=category, color=category))+
  geom_path(linewidth=0.5)+
  geom_point(size=2)+
  scale_color_manual(values=c("grey30", "deepskyblue3"))+
  labs(x="Diet specialization")+
  scale_x_reverse()+
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=0.25)+
  #  annotate("text", x=0.7, y=0.83, label="Generalist", angle=90)+
  #  annotate("text", x=0.7, y=0.92, label="Specialist", angle=90)+
  #$ annotate("segment", x = 0.7, y = 0.89, xend = 0.7, yend = 0.88, size=0.5,
  #        arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
  #annotate("segment", x = 0.7, y = 0.93, xend = 0.7, yend = 0.94, size=0.5,
  #          arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
  theme_classic()+
  theme(axis.title.y=element_blank(), legend.position="none")
# annotations not showing up for some reason rip
diet_point

diet.plot <- ggarrange(diet_bar, diet_point, ncol=1)
ggsave(diet.plot, file="diet.spec.results.png", height=4, width=8)


###### Look at relationship between habitat and diet breadth

# merge habitat and diet breadth for each species
unique.species.diet <- sp_diet %>% group_by(SCIENTIFIC.NAME) %>% summarise(diet.breadth=mean(gini.index)) 
unique.species.habitat <- sp_habitat %>% group_by(SCIENTIFIC.NAME) %>% summarise(habitat.breadth=mean(Habitat_breadth_IUCN))
habitat.diet <- merge(unique.species.diet, unique.species.habitat, by="SCIENTIFIC.NAME") #6724

ggplot(habitat.diet, aes(x=diet.breadth, y=habitat.breadth))+
  geom_point()+
  geom_smooth(method="lm")

cor(habitat.diet$diet.breadth, habitat.diet$habitat.breadth) # negative because they are opposite
# diet and habitat breadth do not correlate very strongly, but they are both showing a signal with urbanization




#################### Figure out how to use taxise to merge data

#res<-taxize::get_gbifid_(global_uniquesp$SCIENTIFIC.NAME, method="backbone") #finds GBIF info for each species 
#all.names<-as.data.frame(matrix(data=NA,nrow=nrow(species.list),ncol=2))
#names(all.names)=c("IUCN_Name","GBIF_Name")
#for (i in 347:length(res)){
#  all.names[i,1]=names(res)[i]
#  
#  if (length(which(res[[i]]$status=="ACCEPTED" & res[[i]]$matchtype=="EXACT"))>0){
#    all.names[i,2]=res[[i]]$species[which(res[[i]]$status=="ACCEPTED" & res[[i]]$matchtype=="EXACT")]
#  }
#  
#  if (length(which(res[[i]]$status=="ACCEPTED" & res[[i]]$matchtype=="EXACT"))==0){
#    all.names[i,2]=res[[i]]$species[which(res[[i]]$status=="SYNONYM" & res[[i]]$matchtype=="EXACT")]
#  }
#  
#  else(next)
#}




######################################
##### Look at proportional richness loss in latitude bins with full data
global_uniquesp <- read.table("global_unique_species.txt", header=TRUE)
length(unique(global_uniquesp$SCIENTIFIC.NAME))
global_uniquesp$abslat <- abs(global_uniquesp$lat)

# Bin by latitude
global_uniquesp <- global_uniquesp %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 90), 
                                                             labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))

total_zone <- global_uniquesp %>% group_by(zone_bin, SCIENTIFIC.NAME) %>% count()

# number of species in each zone
total_zone %>% group_by(zone_bin) %>% count()
# 9106 tripical, 5204 subtropical, 3041 temperate, 376 arctic

total_zones <- global_uniquesp %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

total_zones <- total_zones %>% replace(is.na(.), 0)

total_zones$category <- NA

# label by urban and not urban
for (i in 1:nrow(total_zones)){
  if (total_zones$natural[i] >= 0 & total_zones$urban[i] > 0) {
    total_zones$category[i] <- "both"
  }
  #  else if (total_zones$natural[i] == 0 & total_zones$urban[i] > 0) {
  #    total_zones$category[i] <- "urban.only"
  #  }
  else if (total_zones$natural[i] > 0 & total_zones$urban[i] == 0) {
    total_zones$category[i] <- "not.urban"
  }
}

zzz <- total_zones %>% group_by(zone_bin, category) %>% count() %>% pivot_wider(names_from="category", values_from="n") 

#zzz5 <- total_zones %>% group_by(category) %>% count() %>% pivot_wider(names_from="category", values_from="n") %>% 
# mutate(total=sum(urban, urban.only, not.urban), fraction=urban.only/total)
richness_category <- total_zones %>% group_by(zone_bin, category) %>% count()

## Make plot with overall results of species being lost (because some species lost when merged with habitat or diet data)
total_bar <- ggplot(richness_category, aes(fill=reorder(category, n), y=n, x=reorder(zone_bin, -n))) + 
  scale_fill_manual(labels=c('In natural only', 'In urban'), values=c("deepskyblue3", "grey30"))+
  coord_flip()+
  labs(y="Number of Species")+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(), axis.title.y=element_blank(), legend.title=element_blank(), 
        legend.position = c(0.85, 0.9), legend.text = element_text(size=13), axis.title.x=element_text(size=12),
        axis.text=element_text(size=10), legend.spacing.y = unit(1, 'cm'))+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))
total_bar

# this is the total number of species (not only the ones that matched up)
#ggsave(total_bar, file="total_species_zones.png", height=7, width=4)




######## Put all plots together
library(patchwork)

#habitat_point2 <- emmeans.df.habitat %>% filter(!category=="urban.only") %>% 
#  ggplot(aes(x=zone_bin, y=emmean, group=category, color=category))+
#  geom_point(size=2, position=position_dodge(width=0.2))+
#  geom_line(size=0.5, position=position_dodge(width=0.2))+
#  scale_color_manual(values=c("grey30", "deepskyblue3"))+
#  scale_y_reverse()+
#  labs(y="Log habitat breadth")+
#  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.15, position=position_dodge(width=0.2))+
#  annotate("text", x=0.7, y=2.3, label="Generalist", angle=90)+
#  annotate("text", x=0.7, y=1.55, label="Specialist", angle=90)+
#  annotate("segment", x = 0.7, y = 2.75, xend = 0.7, yend = 2.9, size=0.6,
#           arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
#  annotate("segment", x = 0.7, y = 1, xend = 0.7, yend = 1.15, size=0.5,
#           arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="first"))+
#  coord_cartesian(clip = "off")+
#  theme_classic()+
#  theme(axis.title.x=element_blank(), legend.position="none", axis.title.y=element_text(size=12),
#        axis.text=element_text(size=10))
#habitat_point2
#
#
#diet_point2 <- emmeans.df.diet %>% 
#  ggplot(aes(x=zone_bin, y=emmean, group=category, color=category))+
#  geom_point(size=2, position=position_dodge(width=0.2))+
#  geom_line(linewidth=0.5, position=position_dodge(width=0.2))+
#  scale_color_manual(values=c("grey30", "deepskyblue3"))+
#  labs(y="Diet specialization")+
#  # scale_y_reverse()+
#  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.25, position=position_dodge(width=0.2))+
# # annotate("text", x=0.7, y=0.9, label="Generalist", angle=90)+
#  #annotate("text", x=0.7, y=0.92, label="Specialist", angle=90)+
#  #annotate("segment", x = 0.7, y = 0.89, xend = 0.7, yend = 0.88, size=0.5,
#   #        arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
#  #annotate("segment", x = 0.7, y = 0.93, xend = 0.7, yend = 0.94, size=0.5,
#   #        arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
#  theme_classic()+
#  theme(axis.title.x=element_blank(), legend.position="none", axis.title.y=element_text(size=12),
#        axis.text=element_text(size=10))
## annotations not showing up for some reason rip
#diet_point2

composite_plot <- total_bar | (habitat_point / diet_point) + plot_annotation(tag_levels = "A") + plot_layout(widths=c(2, 1.5))
ggsave(composite_plot, file="full_specialist_results.png", height=9, width=11)


####### Euler plots
# another category with urban only
total_zones$category2 <- NA
for (i in 1:nrow(total_zones)){
  if (total_zones$natural[i] > 0 & total_zones$urban[i] > 0) {
    total_zones$category[i] <- "both"
  }
  else if (total_zones$natural[i] == 0 & total_zones$urban[i] > 0) {
    total_zones$category[i] <- "urban.only"
  }
  else if (total_zones$natural[i] > 0 & total_zones$urban[i] == 0) {
    total_zones$category[i] <- "not.urban"
  }
}


total_zones %>% group_by(zone_bin) %>% count()
richness_category <- total_zones %>% group_by(zone_bin, category) %>% count()
# put it into a dataframe that eulerr will understand
library(eulerr)
euler_trop_df <- c(Ntrop=3556, Utrop=217, "Ntrop&Utrop"=5180, Nsub=1315, Usub=178, 
                   "Nsub&Usub"=3734, Ntemp=608, Utemp=173, "Ntemp&Utemp"=2139, Npol=350, Upol=88, "Npol&Upol"=927) # tropical
euler_subtrop_df <- c(Nsub=1315, Usub=178, "Nsub&Usub"=3734) # subtropical
euler_temp_df <- c(Ntemp=608, Utemp=173, "Ntemp&Utemp"=2139) # temperate
euler_subpol <- c(Npol=350, Upol=88, "Npol&Upol"=927) # subpolar

euler_trop <- euler(euler_trop_df, shape="ellipse")
euler_plot <- plot(euler_trop, fills=c("deepskyblue3", "grey30", "deepskyblue3", "grey30", "deepskyblue3", "grey30", "deepskyblue3", "grey30"), labels="")
png("euler_plot.png")
euler_plot
dev.off()




#######################################

## Look at distribution of different diet guilds
test <- diet_zones %>% pivot_longer(cols=c("natural", "urban"))

ggplot(test)+
  geom_col(aes(x=name, y=value, fill=Diet.5Cat)) +
  facet_wrap(~zone_bin)
# not really seeing anything here



















