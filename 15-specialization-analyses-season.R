############### Specialization analyses for seasonal data ##################
library(tidyverse)
library(terra)
library(taxize)
library(sf)
library(terra)
library(emmeans)

##### Winter #################

winter.model.data <- read.csv("season_modeling_data.csv") %>% filter(season=="winter")

unique(winter.model.data$square)

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c2", "r4c4") # take out r4c1 and r4c4 because they are not in the final modelling data

GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")

datalist.years <- list()
datalist.names <- list()

for (i in 1:length(names)){ # come back to r1c4 (4) - no data in r1c4
  for (j in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% winter.model.data$cell) # filter for cells that are in the final dataset
    dat.filt$SCIENTIFIC.NAME<- as.character(dat.filt$SCIENTIFIC.NAME)
    dat.filt$OBSERVATION.DATE<- as.character(dat.filt$OBSERVATION.DATE)
    dat.filt$OBSERVER.ID<- as.character(dat.filt$OBSERVER.ID)
    dat.filt$SAMPLING.EVENT.IDENTIFIER <- as.character(dat.filt$SAMPLING.EVENT.IDENTIFIER)
    dat.filt$OBSERVATION.COUNT <- as.character(dat.filt$OBSERVATION.COUNT)
    dat.filt$GROUP.IDENTIFIER <- as.character(dat.filt$GROUP.IDENTIFIER)
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  wint_dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  wint_dat_uniquesp$x <- xFromCell(GHSL, wint_dat_uniquesp$cell) # extract the coordinates from the cells
  wint_dat_uniquesp$y <- yFromCell(GHSL, wint_dat_uniquesp$cell)
  wint_dat_uniquesp$square <- names[i]
  datalist.names[[i]] <- wint_dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(wint_dat_uniquesp)
}

winter_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(winter_unique_sp$SCIENTIFIC.NAME)) # 8,724 species
write.table(winter_unique_sp, "winter_unique_species.txt", row.names=FALSE)





############# Summer ##################
summer.model.data <- read.csv("season_modeling_data.csv") %>% filter(season=="summer")


datalist.years <- list()
datalist.names <- list()

for (i in 1:length(names)){ 
  for (j in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), 
                      header=TRUE)
    dat.filt <- dat %>% filter(cell %in% summer.model.data$cell) # filter for cells that are in the final dataset
    dat.filt$SCIENTIFIC.NAME<- as.character(dat.filt$SCIENTIFIC.NAME)
    dat.filt$OBSERVATION.DATE<- as.character(dat.filt$OBSERVATION.DATE)
    dat.filt$OBSERVER.ID<- as.character(dat.filt$OBSERVER.ID)
    dat.filt$SAMPLING.EVENT.IDENTIFIER <- as.character(dat.filt$SAMPLING.EVENT.IDENTIFIER)
    dat.filt$OBSERVATION.COUNT <- as.character(dat.filt$OBSERVATION.COUNT)
    dat.filt$GROUP.IDENTIFIER <- as.character(dat.filt$GROUP.IDENTIFIER)
    datalist.years[[j]] <- dat.filt
  }
  dat2 <- dplyr::bind_rows(datalist.years) # put all years together
  sum_dat_uniquesp <- dat2 %>% 
    distinct(cell, SCIENTIFIC.NAME)
  sum_dat_uniquesp$x <- xFromCell(GHSL, sum_dat_uniquesp$cell) # extract the coordinates from the cells
  sum_dat_uniquesp$y <- yFromCell(GHSL, sum_dat_uniquesp$cell)
  sum_dat_uniquesp$square=names[i]
  datalist.names[[i]] <- sum_dat_uniquesp
  print(paste("finished", names[i]))
  rm(dat)
  rm(dat.filt)
  rm(sum_dat_uniquesp)
}

summer_unique_sp <- dplyr::bind_rows(datalist.names) # put all sections together
length(unique(summer_unique_sp$SCIENTIFIC.NAME)) # 10,723 species
write.table(summer_unique_sp, "summer_unique_species.txt", row.names=FALSE)

##########################

# Merge with trait data
# Combine summer and winter dataframes
winter_uniquesp <- read.table("winter_unique_species.txt", header=TRUE)
winter_uniquesp$season <- "winter"
summer_uniquesp <- read.table("summer_unique_species.txt", header=TRUE)
summer_uniquesp$season <- "summer"
# merge

season_uniquesp <- rbind(winter_uniquesp, summer_uniquesp)

GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif")

length(unique(season_uniquesp$SCIENTIFIC.NAME)) #9373 species


# add in lat long points to more easily bin by latitude
dat_latlong <- st_as_sf(season_uniquesp, coords=c("x", "y"), crs=st_crs(GHSL))
dat_latlong <- st_transform(dat_latlong, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- as.data.frame(dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                                   lat = sf::st_coordinates(.)[,2]))

# bind this with data in other crs
season_uniquesp <- cbind(season_uniquesp, latlong_df[,6:7])

#### Extract urban scores
season_uniquesp$urban <- as.data.frame(terra::extract(GHSL, season_uniquesp[,c(3:4)]))$SMOD_global
test <- season_uniquesp %>% na.omit(urban) # there are no NAs in urban, this is good
# turn urban into 3 categories
season_uniquesp <- season_uniquesp %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), "natural", ifelse(urban==30, "urban", "suburban")))



########## Merge with trait data
habitat <- read.csv("/Volumes/Expansion/eBird/Traits/habitat_breadth.csv")
season_sp_habitat <- merge(season_uniquesp, habitat[,c(4,14)], by.x="SCIENTIFIC.NAME", by.y="Best_guess_binomial")
length(unique(season_sp_habitat$SCIENTIFIC.NAME)) # 8,498 species

season_sp_habitat$abslat <- abs(season_sp_habitat$lat)
# Try with habitat data
season_sp_habitat <- season_sp_habitat %>% mutate(lat_bin = cut(abslat, breaks=abs(c(0, 10, 20, 30, 40, 50, 60, 70, 80))))                                             
# save data with habitat breadth
write.table(season_sp_habitat, "season_habitatbreadth.txt", row.names=F)

######## Diet data 
diet<- read.csv("/Volumes/Expansion/eBird/Traits/EltonTraits/BirdFuncDat_wgini.csv") # load diet data
season_sp_diet <- merge(season_uniquesp, diet[, c(9,42)], by.x="SCIENTIFIC.NAME", by.y="Scientific") # merge with species data
length(unique(season_sp_diet$SCIENTIFIC.NAME)) # 6,904 species

# Bin latitude by 10 degrees
season_sp_diet$abslat <- abs(season_sp_diet$lat)
season_sp_diet <- season_sp_diet %>% mutate(lat_bin = cut(abslat, breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80)))
# save data with diet
write.table(season_sp_diet, "season_dietspec.txt", row.names=F)


############################################




########## Modelling season habitat
season_sp_habitat <- read.table("season_habitatbreadth.txt", header=TRUE) %>% filter(!is.na(Habitat_breadth_IUCN))

# Plot of habitat breadth and urbanization by latitude bin

# Divide into urban only, both, and natural only
#summer_categories <- summer_sp_habitat %>% group_by(lat_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
 # filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  
#summer_categories <- summer_categories %>% replace(is.na(.), 0)

#summer_categories$category <- NA
# to make it simpler I will take out surburban for now
#for (i in 1:nrow(summer_categories)){
 # if (summer_categories$natural[i] > 0 & summer_categories$urban[i] > 0) {
#    summer_categories$category[i] <- "both"
 # }
#  else if (summer_categories$natural[i] > 0 & summer_categories$urban[i] == 0) {
 #   summer_categories$category[i] <- "natural.only"
#  }
 # else if (summer_categories$natural[i] == 0 & summer_categories$urban[i] > 0) {
#    summer_categories$category[i] <- "urban.only"
 # }
#}

# make boxplot
#ggplot(summer_categories)+
 # geom_boxplot(aes(x=lat_bin, y=log(Habitat_breadth_IUCN), fill=category))

#### Try binning by larger categories
season_sp_habitat <- season_sp_habitat %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90), labels=c("Tropical", "Subtropical", "Temperate", "Arctic")))

season_zones <- season_sp_habitat %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN, season) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

season_zones <- season_zones %>% replace(is.na(.), 0)

season_zones$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(season_zones)){
  if (season_zones$natural[i] > 0 & season_zones$urban[i] > 0) {
    season_zones$category[i] <- "In urban"
  }
  else if (season_zones$natural[i] > 0 & season_zones$urban[i] == 0) {
    season_zones$category[i] <- "Not in urban"
  }
  else if (season_zones$natural[i] == 0 & season_zones$urban[i] > 0) {
    season_zones$category[i] <- "urban.only"
  }
}


ggplot(season_zones)+
  geom_boxplot(aes(x=zone_bin, y=log(Habitat_breadth_IUCN), fill=category))

# run an anova
season.habitat.aov <- aov(log(Habitat_breadth_IUCN) ~ zone_bin * category * season, data = season_zones)
summary(season.habitat.aov)
emmeans.results <- emmeans(season.habitat.aov, specs=c("season", "category"), by="zone_bin", facet=TRUE)
plot(emmeans.results)
#### Looks pretty similar to winter and overall

richness_category <- season_zones %>% group_by(zone_bin, category, season) %>% count()
season_bar <-ggplot(richness_category, aes(fill=reorder(category, n), y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('Absent from urban', 'Found in urban'), values=c("deepskyblue3", "black"))+
  geom_bar(position="stack", stat="identity")+
  labs(y="Number of Species")+
  #  coord_flip()+
  theme_classic()+
  facet_wrap(~season)
season_bar

# Plot of habitat breadth means
season.emmeans.df <- as.data.frame(emmeans.results)

ggplot(season.emmeans.df, aes(x=zone_bin, y=emmean, group=category, color=category))+
  geom_point(size=2)+
  geom_line(linewidth=0.5)+
  scale_color_manual(labels=c('In urban', 'Not in urban'), values=c("#000000","#009E73"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.25)+
  facet_wrap(~season)+
  theme_bw()
# specialization measures for urban and not urban are the same for summer and winter
# difference between specialization are definitely decreasing with latitude



##############################################

## Diet
season_sp_diet <- read.table("season_dietspec.txt", header=TRUE) %>% filter(!is.na(gini.index))

# Plot of habitat breadth and urbanization by latitude bin

# Divide into urban only, both, and natural only
#summer_categories <- summer_sp_habitat %>% group_by(lat_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN) %>% count(.drop=FALSE) %>% 
# filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  
#summer_categories <- summer_categories %>% replace(is.na(.), 0)

#summer_categories$category <- NA
# to make it simpler I will take out surburban for now
#for (i in 1:nrow(summer_categories)){
# if (summer_categories$natural[i] > 0 & summer_categories$urban[i] > 0) {
#    summer_categories$category[i] <- "both"
# }
#  else if (summer_categories$natural[i] > 0 & summer_categories$urban[i] == 0) {
#   summer_categories$category[i] <- "natural.only"
#  }
# else if (summer_categories$natural[i] == 0 & summer_categories$urban[i] > 0) {
#    summer_categories$category[i] <- "urban.only"
# }
#}

# make boxplot
#ggplot(summer_categories)+
# geom_boxplot(aes(x=lat_bin, y=log(Habitat_breadth_IUCN), fill=category))

#### Try binning by larger categories
season_sp_diet <- season_sp_diet %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 66.5, 90), labels=c("Tropical", "Subtropical", "Temperate", "Arctic")))

season_zones_diet <- season_sp_diet %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, gini.index, season) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

season_zones_diet <- season_zones_diet %>% replace(is.na(.), 0)

season_zones_diet$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(season_zones_diet)){
  if (season_zones_diet$natural[i] >= 0 & season_zones_diet$urban[i] > 0) {
    season_zones_diet$category[i] <- "In urban"
  }
  else if (season_zones_diet$natural[i] > 0 & season_zones_diet$urban[i] == 0) {
    season_zones_diet$category[i] <- "Not in urban"
  }
  # else if (season_zones_diet$natural[i] == 0 & season_zones_diet$urban[i] > 0) {
  #  season_zones_diet$category[i] <- "urban.only"
  #}
}


ggplot(season_zones_diet)+
  geom_boxplot(aes(x=zone_bin, y=gini.index, fill=category))

# run an anova
season.diet.aov <- aov(gini.index ~ zone_bin * category * season, data = season_zones_diet)
summary(season.diet.aov)
emmeans.results <- emmeans(season.diet.aov, specs=c("season", "category"), by="zone_bin", facet=TRUE)
plot(emmeans.results)
#### Looks pretty similar to winter and overall

richness_category <- season_zones_diet %>% group_by(zone_bin, category, season) %>% count()
season_bar <-ggplot(richness_category, aes(fill=reorder(category, n), y=n, x=zone_bin)) + 
  scale_fill_manual(labels=c('Absent from urban', 'Found in urban'), values=c("deepskyblue3", "black"))+
  geom_bar(position="stack", stat="identity")+
  labs(y="Number of Species")+
  #  coord_flip()+
  theme_classic()+
  facet_wrap(~season)
season_bar

# Plot of diet breadth means
season.emmeans.df <- as.data.frame(emmeans.results)

ggplot(season.emmeans.df, aes(x=zone_bin, y=emmean, group=category, color=category))+
  geom_point(size=2)+
  geom_line(linewidth=0.5)+
  scale_color_manual(labels=c('In urban', 'Not in urban'), values=c("#000000","#009E73"))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.25)+
  facet_wrap(~season)+
  theme_bw()
# this looks wrong because they are negative? but maybe because it
