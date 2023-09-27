########### Script for analyzing whether specialist species are being disproportionately lost at lower latitudes
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)
library(emmeans)
## First I want to thin the data so that for each cell, there is only one row for each species 

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
          "r4c2", "r4c4") # take out r4c1 and r4c4 because they are not in the final modelling data


model.data <- read.csv("modeling_data.csv") 
# also going to load this because I want to filter for cells that are in the final dataset to save space
unique(model.data$square)
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

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
write.table(global_unique_sp, "global_unique_species.txt", row.names=FALSE)


###################################################################################################

GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif")
global_uniquesp <- read.table("global_unique_sp.txt", header=TRUE) # this worked better than csv
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



### Merge with trait data

############ Habitat data
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")
sp_habitat <- merge(global_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
length(unique(sp_habitat$SCIENTIFIC.NAME)) # 8,498 species

sp_habitat$abslat <- abs(sp_habitat$lat)
# Try with habitat data
sp_habitat <- sp_habitat %>% mutate(lat_bin = cut(abslat, breaks=abs(c(0, 10, 20, 30, 40, 50, 60, 70, 80))))                                             
# save data with habitat breadth
write.table(sp_habitat, "unique_sp_habitatbreadth.txt", row.names=F)

######## Diet data 
diet<- read.csv("/Volumes/Expansion/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
sp_diet <- merge(global_uniquesp, diet[, c(9,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific") # merge with species data
length(unique(sp_diet$SCIENTIFIC.NAME)) # 6,904 species

# Bin latitude by 10 degrees
sp_diet$abslat <- abs(sp_diet$lat)
sp_diet <- sp_diet %>% mutate(lat_bin = cut(abslat, breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80)))
# save data with diet
write.table(sp_diet, "unique_sp_dietspec.txt", row.names=F)





#############################################################




## Habitat data 
sp_habitat <- read.table("unique_sp_habitatbreadth.txt", header=T)
# see if there are more specialists at low latitudes
# boxplot of habitat specialization
hist(log(sp_habitat$Habitat_breadth_IUCN)) # definitely looks pretty log normal
hist(sp_habitat$Habitat_breadth_IUCN)

ggplot(sp_habitat)+
  geom_boxplot(aes(x=lat_bin, y=log(Habitat_breadth_IUCN), fill=urban2))
# looks like there could be a pattern here - habitat breadth decreases at 
# low latitude bins in natural but not urban areas
ggplot(sp_habitat)+
  geom_boxplot(aes(x=lat_bin, y=log(Habitat_breadth_IUCN)))

habitat.aov <- aov(Habitat_breadth_IUCN ~ lat_bin, data = sp_habitat)
summary(habitat.aov)
# there does seem to be a significant difference
library(emmeans)
emmeans.results1 <- emmeans(habitat.aov, specs="lat_bin")
plot(emmeans.results1)


# anova
habitat.aov <- aov(Habitat_breadth_IUCN ~ lat_bin * urban2, data = sp_habitat)
summary(habitat.aov)
plot(habitat.aov)
# there does seem to be a significant difference
emmeans.results1 <- emmeans(habitat.aov, specs="urban2", by="lat_bin")
plot(emmeans.results1) # wow this is super interesting too! The difference definitely decreases

## Do a big grouping by species and latitude bin and label by both urban and non-urban, just urban, or just non-urban
# to make it simpler I will take out suburban for now
birds_test <- sp_habitat %>% group_by(lat_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

birds_test <- birds_test %>% replace(is.na(.), 0)

birds_test$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(birds_test)){
  if (birds_test$natural[i] >= 0 & birds_test$urban[i] > 0) {
    birds_test$category[i] <- "urban"
  }
  else if (birds_test$natural[i] > 0 & birds_test$urban[i] == 0) {
    birds_test$category[i] <- "natural.only"
}
}

birds_test %>% group_by(category) %>% count()
# there are less species that are natural only

# make boxplot of habitat breadth for each latitude bin with habitat breadth and urbanization
ggplot(birds_test)+
  geom_boxplot(aes(x=lat_bin, y=Habitat_breadth_IUCN, fill=category))

# run anova
habitat.aov2 <- aov(Habitat_breadth_IUCN ~ lat_bin * category, data = birds_test)
summary(habitat.aov2)

emmeans.results2 <- emmeans(habitat.aov2, specs="category", by="lat_bin")
plot(emmeans.results2)
# this is super cool!


##### Bin by larger groupings
sp_habitat <- sp_habitat %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90)))

# run anova on the raw habitat breadth with these larger zones
habitat.aov3 <- aov(Habitat_breadth_IUCN ~ zone_bin * urban2, data = sp_habitat)
summary(habitat.aov3) # significant interaction

emmeans.results3 <- emmeans(habitat.aov3, specs="urban2", by="zone_bin")
plot(emmeans.results3)
# the difference is way larger in the tropics!

birds_zones <- sp_habitat %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

birds_zones <- birds_zones %>% replace(is.na(.), 0)

birds_zones$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(birds_zones)){
  if (birds_zones$natural[i] >= 0 & birds_zones$urban[i] > 0) {
    birds_zones$category[i] <- "urban"
  }
  else if (birds_zones$natural[i] > 0 & birds_zones$urban[i] == 0) {
    birds_zones$category[i] <- "lost.in.urban"
  }
}

habitat.aov4 <- aov(log(Habitat_breadth_IUCN) ~ zone_bin * category, data = birds_zones)
summary(habitat.aov4)

emmeans.results4 <- emmeans(habitat.aov4, specs="category", by="zone_bin")
plot(emmeans.results4)

# Plot of richness of different categories
richness_category <- birds_zones %>% group_by(zone_bin, category) %>% count()
ggplot(richness_category, aes(fill=category, y=n, x=zone_bin)) + 
  geom_bar(position="stack", stat="identity")
# this is exactly what I was trying to show! Didn't know it would be so clear

# Make ridgeline plot
# Jenn asked me to do this but it is not telling that much to be honest
library(ggridges)
ggplot(birds_zones, aes(y = zone_bin, x = Habitat_breadth_IUCN, fill = category)) +
  geom_density_ridges() +
  theme_ridges() 




###########################################





### Diet specialization
sp_diet <- read.table("unique_sp_dietspec.txt", header=T)
# see if there are more specialists at low latitudes
sp_diet %>% group_by(lat_bin) %>% summarise(mean_diet = mean(gini.index))
# boxplot of specialization
ggplot(sp_diet)+
  geom_boxplot(aes(x=lat_bin, y=gini.index, fill=urban2))

# ok they look pretty similar lol
# anova
diet.aov1 <- aov(gini.index ~ lat_bin * urban2, data = sp_diet)
summary(diet.aov1) # interaction is significant
# look at contrasts

emmeans.results1 <- emmeans(diet.aov1, specs="urban2", by="lat_bin")
plot(emmeans.results1)


## Group by whether they are found in urban, non-urban, etc.
birds_diet <- sp_diet %>% group_by(lat_bin, SCIENTIFIC.NAME, urban2, gini.index) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

birds_diet <- birds_diet %>% replace(is.na(.), 0)

birds_diet$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(birds_diet)){
  if (birds_diet$natural[i] >= 0 & birds_diet$urban[i] > 0) {
    birds_diet$category[i] <- "urban"
  }
  else if (birds_diet$natural[i] > 0 & birds_diet$urban[i] == 0) {
    birds_diet$category[i] <- "natural.only"
  }
}

ggplot(birds_diet)+
  geom_boxplot(aes(x=lat_bin, y=gini.index, fill=category))

# run anova
diet.aov2 <- aov(gini.index ~ lat_bin * category, data = birds_diet)
summary(diet.aov2) # interaction is not significant
# look at contrasts

emmeans.results2 <- emmeans(diet.aov2, specs="category", by="lat_bin")
plot(emmeans.results2)


####### Bin by larger categories

sp_diet <- sp_diet %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90)))

# run anova on the raw habitat breadth with these larger zones
diet.aov3 <- aov(gini.index ~ zone_bin * urban2, data = sp_diet)
summary(diet.aov3) # significant interaction

emmeans.results3 <- emmeans(diet.aov3, specs="urban2", by="zone_bin")
plot(emmeans.results3)
# the difference is way larger in the tropics!

diet_zones <- sp_diet %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, gini.index) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

diet_zones <- diet_zones %>% replace(is.na(.), 0)

diet_zones$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(diet_zones)){
  if (diet_zones$natural[i] >= 0 & diet_zones$urban[i] > 0) {
    diet_zones$category[i] <- "urban"
  }
  else if (diet_zones$natural[i] > 0 & diet_zones$urban[i] == 0) {
    diet_zones$category[i] <- "natural.only"
  }
}

diet.aov4 <- aov(gini.index ~ zone_bin * category, data = diet_zones)
summary(diet.aov4)

emmeans.results4 <- emmeans(diet.aov4, specs="category", by="zone_bin")
plot(emmeans.results4)





###########################
# Change databases so that generalist and specialist are discrete so I can merge it with other data and plot
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")
hist(habitat$Habitat_breadth_IUCN)


# Let's say everything under 4 is a specialist (this is completely arbitrary)



#################### Figure out how to use taxise to merge data

res<-taxize::get_gbifid_(global_uniquesp$SCIENTIFIC.NAME, method="backbone") #finds GBIF info for each species 
all.names<-as.data.frame(matrix(data=NA,nrow=nrow(species.list),ncol=2))
names(all.names)=c("IUCN_Name","GBIF_Name")
for (i in 347:length(res)){
  all.names[i,1]=names(res)[i]
  
  if (length(which(res[[i]]$status=="ACCEPTED" & res[[i]]$matchtype=="EXACT"))>0){
    all.names[i,2]=res[[i]]$species[which(res[[i]]$status=="ACCEPTED" & res[[i]]$matchtype=="EXACT")]
  }
  
  if (length(which(res[[i]]$status=="ACCEPTED" & res[[i]]$matchtype=="EXACT"))==0){
    all.names[i,2]=res[[i]]$species[which(res[[i]]$status=="SYNONYM" & res[[i]]$matchtype=="EXACT")]
  }
  
  else(next)
}








