########### Script for analyzing whether specialist species are being disproportionately lost at lower latitudes
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)
library(emmeans)
library(ggpubr)
library(grid)
library(cowplot)
library(tidyterra)
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

write.table(global_uniquesp, "global_unique_species.txt", row.names=FALSE)




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


## Habitat data, filter out NAs for habitat
sp_habitat <- read.table("unique_sp_habitatbreadth.txt", header=T) %>% filter(!is.na(Habitat_breadth_IUCN))
length(unique(sp_habitat$SCIENTIFIC.NAME))
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

################ Analyses with the 10 degree habitat bins
# anova
#habitat.aov <- aov(Habitat_breadth_IUCN ~ lat_bin * urban2, data = sp_habitat)
#summary(habitat.aov)
##plot(habitat.aov)
## there does seem to be a significant difference
#emmeans.habitat1 <- emmeans(habitat.aov, specs="urban2", by="lat_bin")
#plot(emmeans.habitat1) # wow this is super interesting too! The difference definitely decreases
#
### Do a big grouping by species and latitude bin and label by both urban and non-urban, just urban, or just non-urban
## to make it simpler I will take out suburban for now
#birds_test <- sp_habitat %>% group_by(lat_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
##  filter(!urban2=="suburban") %>% 
#  pivot_wider(names_from="urban2", values_from="n")  
#
#birds_test <- birds_test %>% replace(is.na(.), 0)
#
#birds_test$category <- NA
## to make it simpler I will take out surburban for now
#for (i in 1:nrow(birds_test)){
#  if (birds_test$natural[i] >= 0 & birds_test$suburban[i] >= 0 & birds_test$urban[i] > 0) {
#    birds_test$category[i] <- "in.urban"
#  }
#  else if (birds_test$natural[i] >= 0 & birds_test$suburban[i] > 0 & birds_test$urban[i] == 0){
#    birds_test$category[i] <- "in.suburban"
#  }
#  else if (birds_test$natural[i] > 0 & birds_test$suburban[i] == 0 & birds_test$urban[i] == 0) {
#    birds_test$category[i] <- "natural.only"
#  }
#  
#}
#
#birds_test %>% group_by(category) %>% count()
## there are less species that are natural only
#
## make boxplot of habitat breadth for each latitude bin with habitat breadth and urbanization
#ggplot(birds_test)+
#  geom_boxplot(aes(x=lat_bin, y=Habitat_breadth_IUCN, fill=category))
#
## run anova
#habitat.aov2 <- aov(Habitat_breadth_IUCN ~ lat_bin * category, data = birds_test)
#summary(habitat.aov2)
#
#emmeans.habitat2 <- emmeans(habitat.aov2, specs="category", by="lat_bin")
#plot(emmeans.habitat2)
# this is super cool!


######### Bin by larger groupings
sp_habitat <- sp_habitat %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90), labels=c("Tropical", "Subtropical", "Temperate", "Arctic")))

# run anova on the raw habitat breadth with these larger zones
habitat.aov3 <- aov(Habitat_breadth_IUCN ~ zone_bin * urban2, data = sp_habitat)
summary(habitat.aov3) # significant interaction

# look at just difference in specialization with latitude
emmeans(habitat.aov3, specs="zone_bin")

emmeans.habitat3 <- emmeans(habitat.aov3, specs="urban2", by="zone_bin")
plot(emmeans.habitat3)
# the difference is way larger in the tropics!



birds_zones <- sp_habitat %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
 # filter(!urban2=="suburban") %>% 
  pivot_wider(names_from="urban2", values_from="n")  

birds_zones <- birds_zones %>% replace(is.na(.), 0)

birds_zones$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(birds_zones)){
  if (birds_zones$natural[i] >= 0 & birds_zones$suburban[i] >= 0 & birds_zones$urban[i] > 0) {
    birds_zones$category[i] <- "in.urban"
  }
  else if (birds_zones$natural[i] >= 0 & birds_zones$suburban[i] > 0 & birds_zones$urban[i] == 0){
    birds_zones$category[i] <- "in.suburban"
  }
  else if (birds_zones$natural[i] > 0 & birds_zones$suburban[i] == 0 & birds_zones$urban[i] == 0) {
    birds_zones$category[i] <- "natural.only"
  }
  
}



habitat.aov4 <- aov(log(Habitat_breadth_IUCN) ~ zone_bin * category, data = birds_zones)
summary(habitat.aov4)

emmeans.habitat4 <- emmeans(habitat.aov4, specs="category", by="zone_bin")
emmeans(habitat.aov4, pairwise~zone_bin, by="category") # all means are significantly different from one another across zones

emmeans.habitat4
plot(emmeans.habitat4) # now that I fixed that error there is a significant interaction
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
  mutate(category = factor(category, levels = c('natural.only', 'in.suburban', 'in.urban'), ordered = TRUE)) %>%
  ggplot(aes(fill=category, y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('Found in natural only', 'Found in suburban', 'Found in urban'), values=c("deepskyblue3", "purple1", "black"))+
  labs(y="Number of Species")+
#  coord_flip()+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title=element_blank(),legend.position = c(.95, .75),
        legend.justification = c("right", "bottom"),
       legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
habitat_bar
#  theme(axis.title.x=element_blank(), legend.position="none")

### Plot emmeans using ggplot
emmeans.df.habitat <- as.data.frame(emmeans.habitat4)


habitat_point <-ggplot(emmeans.df.habitat, aes(x=zone_bin, y=emmean, group=category, color=category))+
  geom_point(size=2, position=position_dodge(width=0.2))+
  geom_line(size=0.5, position=position_dodge(width=0.2))+
  scale_color_manual(values=c("purple1", "deepskyblue3", "black"))+
  labs(y="Log habitat breadth")+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.15, position=position_dodge(width=0.2))+
  annotate("text", x=0.7, y=2.5, label="Generalist", angle=90)+
  annotate("text", x=0.7, y=1.4, label="Specialist", angle=90)+
  annotate("segment", x = 0.7, y = 2.75, xend = 0.7, yend = 2.9, size=0.6,
           arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
  annotate("segment", x = 0.7, y = 1, xend = 0.7, yend = 1.15, size=0.5,
           arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="first"))+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="none")
habitat_point
 # theme(axis.title.x=element_blank(), legend.title=element_blank(),legend.position = c(.95, .1),
  #      legend.justification = c("right", "bottom"),
   #     legend.box.just = "right",
    #    legend.margin = margin(6, 6, 6, 6))

habitat_plot <- ggarrange(habitat_bar, habitat_point, ncol=1)
habitat_plot
ggsave(habitat_plot, file="pecialistHabitatResults.png")
# Try it as an inset
#ggdraw(habitat_point)+
#  draw_plot({habitat_bar}, width=0.4, height=0.4, x=0.58, y=0.05)


# take a closer look at the species that are lost from urban areas
lost <- birds_zones %>% filter(category=="natural.only") # 2665 of 13,764
sub.only <- birds_zones %>% filter(category=="in.suburban") # 1977
in.urb <- birds_zones %>% filter(category=="in.urban")

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

length(unique(sp_diet$SCIENTIFIC.NAME))
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

emmeans.diet1 <- emmeans(diet.aov1, specs="urban2", by="lat_bin")
plot(emmeans.diet1)


## Group by 10 degree latitude bins
#birds_diet <- sp_diet %>% group_by(lat_bin, SCIENTIFIC.NAME, urban2, gini.index) %>% count(.drop=FALSE) %>% 
#  #filter(!urban2=="suburban") %>% 
#  pivot_wider(names_from="urban2", values_from="n")  
#
#birds_diet <- birds_diet %>% replace(is.na(.), 0)
#
#birds_diet$category <- NA
## to make it simpler I will take out surburban for now
#for (i in 1:nrow(birds_diet)){
#  if (birds_diet$natural[i] >= 0 & birds_diet$suburban[i] >= 0 & birds_diet$urban[i] > 0) {
#    birds_diet$category[i] <- "in.urban"
#  }
#  else if (birds_diet$natural[i] >= 0 & birds_diet$suburban[i] > 0 & birds_diet$urban[i] == 0) {
#    birds_diet$category[i] <- "in.suburban"
#  }
#  
#  else if (birds_diet$natural[i] > 0 & birds_diet$suburban[i] == 0 & birds_diet$urban[i] == 0) {
#    birds_diet$category[i] <- "natural.only"
#  }
#}
#
#
#ggplot(birds_diet)+
#  geom_boxplot(aes(x=lat_bin, y=gini.index, fill=category))
#
## run anova
#diet.aov2 <- aov(gini.index ~ lat_bin * category, data = birds_diet)
#summary(diet.aov2) # interaction is not significant
## look at contrasts
#
#emmeans.diet2 <- emmeans(diet.aov2, specs="category", by="lat_bin")
#plot(emmeans.diet2)
#

####### Bin by larger categories

sp_diet <- sp_diet %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90), labels=c("Tropical", "Subtropical", "Temperate", "Arctic")))


# run anova on the raw habitat breadth with these larger zones
diet.aov3 <- aov(gini.index ~ zone_bin * urban2, data = sp_diet)
summary(diet.aov3) # significant interaction

emmeans.diet3 <- emmeans(diet.aov3, specs="urban2", by="zone_bin")
plot(emmeans.diet3)
# the difference is way larger in the tropics!

diet_zones <- sp_diet %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, gini.index) %>% count(.drop=FALSE) %>% 
#  filter(!urban2=="suburban") %>%
  pivot_wider(names_from="urban2", values_from="n")  

diet_zones <- diet_zones %>% replace(is.na(.), 0)

diet_zones$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(diet_zones)){
  if (diet_zones$natural[i] >= 0 & diet_zones$suburban[i] >= 0 & diet_zones$urban[i] > 0) {
    diet_zones$category[i] <- "in.urban"
  }
  else if (diet_zones$natural[i] >= 0 & diet_zones$suburban[i] > 0 & diet_zones$urban[i] == 0) {
    diet_zones$category[i] <- "in.suburban"
  }
  
  else if (diet_zones$natural[i] > 0 & diet_zones$suburban[i] == 0 & diet_zones$urban[i] == 0) {
    diet_zones$category[i] <- "natural.only"
  }
}


unique(diet_zones$zone_bin)
diet.aov4 <- aov(gini.index ~ zone_bin * category, data = diet_zones)
summary(diet.aov4)

emmeans.diet4 <- emmeans(diet.aov4, specs="category", by="zone_bin")
emmeans(diet.aov4, pairwise~zone_bin, by="category")
plot(emmeans.diet4)

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
  mutate(category = factor(category, levels = c('natural.only', 'in.suburban', 'in.urban'), ordered = TRUE)) %>%
  ggplot(aes(fill=category, y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('Found in natural only', 'Found in suburban', 'Found in urban'), values=c("deepskyblue3", "purple1", "black"))+
  geom_bar(position="stack", stat="identity")+
  labs(y="Number of Species")+
#  coord_flip()+
  theme_bw()+
#  theme(axis.title.x=element_blank(), legend.position="none")
  theme(axis.title.x=element_blank(), legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(.95, .75),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
diet_bar




emmeans.df.diet <- as.data.frame(emmeans.diet4)
diet_point <- ggplot(emmeans.df.diet, aes(x=zone_bin, y=emmean, group=category, color=category))+
  geom_point(size=2, position=position_dodge(width=0.2))+
  geom_line(linewidth=0.5, position=position_dodge(width=0.2))+
  scale_color_manual(values=c("purple1","black","deepskyblue3"))+
  labs(y="Diet specialization")+
  scale_y_reverse()+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.25, position=position_dodge(width=0.2))+
  annotate("text", x=0.7, y=0.86, label="Generalist", angle=90)+
  annotate("text", x=0.7, y=0.915, label="Specialist", angle=90)+
  annotate("segment", x = 0.7, y = 0.835, xend = 0.7, yend = 0.845, size=0.5,
           arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="first"))+
  annotate("segment", x = 0.7, y = 0.93, xend = 0.7, yend = 0.94, size=0.5,
           arrow = arrow(type = "open", length = unit(0.05, "npc"), ends="last"))+
  theme_classic()+
  theme(axis.title.x=element_blank(), legend.title=element_blank(), legend.position="none")
# annotations not showing up for some reason rip
diet_point

ggarrange(habitat_bar, diet_bar, habitat_point, diet_point)



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




######################################
##### Look at proportional richness loss in latitude bins with full data
global_uniquesp <- read.table("global_unique_species.txt", header=TRUE)
length(unique(global_uniquesp$SCIENTIFIC.NAME))
global_uniquesp$abslat <- abs(global_uniquesp$lat)

# Bin by latitude
global_uniquesp <- global_uniquesp %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90), 
                                             labels=c("Tropical", "Subtropical", "Temperate", "Arctic")))

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
    total_zones$category[i] <- "urban"
  }
  else if (total_zones$natural[i] > 0 & total_zones$urban[i] == 0) {
    total_zones$category[i] <- "not.urban"
  }
}

zzz <- total_zones %>% group_by(zone_bin, category) %>% count() %>% pivot_wider(names_from="category", values_from="n") %>% 
  mutate(total=sum(urban, not.urban), fraction=not.urban/total)
richness_category <- total_zones %>% group_by(zone_bin, category) %>% count()
## Make plot with overall results of species being lost (because some species lost when merged with habitat or diet data)
total_bar <- ggplot(richness_category, aes(fill=reorder(category, n), y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('Absent from urban', 'Found in urban'), values=c("deepskyblue3", "black"))+
  labs(y="Number of Species")+
  #  coord_flip()+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title=element_blank(),legend.position = c(.95, .75),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
total_bar

##### This will probably be the final figure?
ggarrange(total_bar, habitat_point, diet_point, nrow=3, align="hv")





