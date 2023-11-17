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
#write.table(winter_unique_sp, "winter_unique_species.txt", row.names=FALSE)





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
#write.table(summer_unique_sp, "summer_unique_species.txt", row.names=FALSE)

##########################

# Merge with trait data
# Combine summer and winter dataframes
#winter_uniquesp <- read.table("winter_unique_species.txt", header=TRUE)
winter_uniquesp$season <- "winter"
#summer_uniquesp <- read.table("summer_unique_species.txt", header=TRUE)
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
season_uniquesp$abslat <- abs(season_uniquesp$lat)
#### Extract urban scores
season_uniquesp$urban <- as.data.frame(terra::extract(GHSL, season_uniquesp[,c(3:4)]))$SMOD_global
test <- season_uniquesp %>% na.omit(urban) # there are no NAs in urban, this is good
# turn urban into 3 categories
season_uniquesp <- season_uniquesp %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), "natural", ifelse(urban==30, "urban", "suburban")))

write.table(season_uniquesp, "season_unique_species.txt", row.names=FALSE)



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
season_sp_habitat <- season_sp_habitat %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 90), labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))

season_zones <- season_sp_habitat %>% group_by(zone_bin, SCIENTIFIC.NAME, urban2, Habitat_breadth_IUCN, season) %>% count(.drop=FALSE) %>% 
  filter(!urban2=="suburban") %>% pivot_wider(names_from="urban2", values_from="n")  

season_zones <- season_zones %>% replace(is.na(.), 0)

season_zones$category <- NA
# to make it simpler I will take out surburban for now
for (i in 1:nrow(season_zones)){
  if (season_zones$natural[i] >= 0 & season_zones$urban[i] > 0) {
    season_zones$category[i] <- "In urban"
  }
 # else if (season_zones$natural[i] > 0 & season_zones$urban[i] == 0) {
  #  season_zones$category[i] <- "Not in urban"
  #}
  else if (season_zones$natural[i] > 0 & season_zones$urban[i] == 0) {
    season_zones$category[i] <- "Not in urban"
  }
}


ggplot(season_zones)+
  geom_boxplot(aes(x=zone_bin, y=log(Habitat_breadth_IUCN), fill=category))
library(ggeffects)
# run an anova
season.habitat.aov <- aov(log(Habitat_breadth_IUCN) ~ zone_bin * category * season, data = season_zones)
summary(season.habitat.aov)
emmeans.results <- emmeans(season.habitat.aov, specs=c("season", "category"), by="zone_bin", facet=TRUE)
hypothesis_test(season.habitat.aov, terms=c("season", "category", "zone_bin"))
?hypothesis_test
plot(emmeans.results)
#### Looks pretty similar to winter and overall
# the average values of specialization are the same in winter and summer in each latitude bin


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

season.habitat.plot <- ggplot(season.emmeans.df, mapping=aes(x=zone_bin, y=emmean, group=interaction(category, season), color=category, shape=season))+
  geom_point(size=2)+
  scale_y_reverse()+
  geom_line(linewidth=0.5, aes(group=interaction(category, season), linetype=season))+
  scale_color_manual(labels=c('In urban', 'Not in urban'), values=c("#000000","deepskyblue3"))+
  scale_linetype_manual(labels=c("Summer", "Winter"), values=c(1,2))+
  scale_shape_manual(labels=c("Summer", "Winter"), values=c(15,17))+
  labs(y="Habitat breadth")+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.25)+
  theme_classic()+
  theme(axis.title.x=element_blank(), legend.title=element_blank(), legend.position=c(0.8, 0.8))
# specialization measures for urban and not urban are the same for summer and winter
# difference between specialization are definitely decreasing with latitude
ggsave(season.habitat.plot, file="season.habitat.plot.png", height=6, width=9)

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
season_sp_diet <- season_sp_diet %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 90), labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))

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
hypothesis_test(season.diet.aov, terms=c("season", "category", "zone_bin"))
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




season.diet.plot <- ggplot(season.emmeans.df, mapping=aes(x=zone_bin, y=emmean, group=interaction(category, season), color=category, shape=season))+
  geom_point(size=2)+
 # scale_y_reverse()+
  geom_line(linewidth=0.5, aes(group=interaction(category, season), linetype=season))+
  scale_color_manual(labels=c('In urban', 'Not in urban'), values=c("#000000","deepskyblue3"))+
  scale_linetype_manual(labels=c("Summer", "Winter"), values=c(1,2))+
  scale_shape_manual(labels=c("Summer", "Winter"), values=c(15,17))+
  labs(y="Diet breadth")+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.25)+
  theme_classic()+
  theme(axis.title.x=element_blank(), legend.title=element_blank(), legend.position="none")

season.diet.plot
ggsave(season.diet.plot, file="diet.specialization.season.png", height=5, width=9)






#################### Plot full species list (not just merged) #######################
season_uniquesp <- read.table("season_unique_species.txt", header=TRUE)

season_uniquesp <- season_uniquesp %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 90), 
                                                             labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))


total_zone <- season_uniquesp %>% group_by(zone_bin, season, SCIENTIFIC.NAME) %>% count()

# number of species in each zone
total_zone %>% group_by(zone_bin) %>% count()

total_zones <- season_uniquesp %>% group_by(zone_bin, SCIENTIFIC.NAME, season, urban2) %>% count(.drop=FALSE) %>% 
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

# seeing how many species switched from not urban to urban in summer to winter
wide_zones <- total_zones[, -c(4,5)] %>% pivot_wider(values_from="category", names_from="season") %>% drop_na() %>% filter(!summer==winter) %>% 
  filter(zone_bin=="Subpolar", summer=="both")






zzz <- total_zones %>% group_by(zone_bin, category) %>% count() %>% pivot_wider(names_from="category", values_from="n") 

#zzz5 <- total_zones %>% group_by(category) %>% count() %>% pivot_wider(names_from="category", values_from="n") %>% 
# mutate(total=sum(urban, urban.only, not.urban), fraction=urban.only/total)
richness_category <- total_zones %>% group_by(zone_bin, season, category) %>% count()

## Make plot with overall results of species being lost (because some species lost when merged with habitat or diet data)
total_bar <- ggplot(richness_category, aes(fill=reorder(category, n), y=n, x=reorder(zone_bin, -n))) + 
  scale_fill_manual(labels=c('In natural only', 'In urban'), values=c("deepskyblue3", "grey30"))+
  #  coord_flip()+
  labs(y="Number of Species")+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(), 
        legend.position = c(0.85, 0.8), legend.text = element_text(size=13), axis.title.y=element_text(size=12),
        axis.text=element_text(size=10), legend.spacing.y = unit(1, 'cm'))+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))+
  facet_wrap(~season)
total_bar
ggsave(total_bar, file="season_geographiczone.png", height=6, width=10)
# this is the total number of species (not only the ones that matched up)
proportion <- richness_category %>% pivot_wider(names_from="category", values_from="n") %>% mutate(total=sum(both+not.urban), prop = not.urban/total)
# although there are less species overall in summer, the proportional richness loss is the same, which is likely why specialization is the same




################# Combine data with migration data to see if migratory birds are leaving low latitudes in summer
avonet <- read.csv("/Users/jorygriffith/Desktop/Avonet/ELEData/TraitData/AVONET2_eBird.csv")
# for migration 1=sedentary, 2 = partially migratory, and 3 = migratory
season_uniquesp <- read.table("season_unique_species.txt", header=TRUE)
# merge with specialization data
season_migration <- merge(season_uniquesp, avonet[,c(1,28)], by.x="SCIENTIFIC.NAME", by.y="Species2")
length(unique(season_migration$SCIENTIFIC.NAME)) # 7992 unique species

season_migration <- season_migration %>% mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 90), 
                                                             labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))


total_zones <- season_migration %>% group_by(zone_bin, season, Migration, SCIENTIFIC.NAME) %>% count()

total_zones2 <- total_zones %>% group_by(zone_bin, season, Migration) %>% count()
total_zones2$Migration <- as.factor(total_zones2$Migration)


## Make plot with overall results of species being lost (because some species lost when merged with habitat or diet data)
total_bar <- ggplot(total_zones2, aes(fill=Migration, y=n, x=zone_bin)) + 
  #  coord_flip()+
  labs(y="Number of Species")+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))+
  facet_wrap(~season)
total_bar
# there aren't any crazy patterns here, although there are more migrants at higher latitudes in the summer

# plot proportion
props <- total_zones2 %>% group_by(zone_bin, season, urban2) %>% mutate(total=sum(n), proportion=n/total)

ggplot(props, aes(x=zone_bin, y=proportion, fill=Migration))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~season)


### Now separate by urbanization level
total_zones <- season_migration %>% group_by(zone_bin, season, Migration, SCIENTIFIC.NAME, urban2) %>% count()

total_zones2 <- total_zones %>% group_by(zone_bin, season, Migration, urban2) %>% count()
total_zones2$Migration <- as.factor(total_zones2$Migration)

total_bar <- ggplot(total_zones2, aes(fill=Migration, y=n, x=zone_bin)) + 
  #  coord_flip()+
  labs(y="Number of Species")+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))+
  facet_grid(urban2~season)
total_bar


#### Now separate by in urban and not in urban

total_zones <- season_migration %>% group_by(zone_bin, SCIENTIFIC.NAME, season, urban2, Migration) %>% count(.drop=FALSE) %>% 
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

richness_category <- total_zones %>% group_by(zone_bin, season, category, Migration) %>% count()

total_bar <- ggplot(richness_category, aes(fill=reorder(category, n), y=n, x=zone_bin)) + 
  #  coord_flip()+
  labs(y="Number of Species")+
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))+
  facet_grid(Migration~season)
total_bar

props <- richness_category %>% group_by(zone_bin, season, category) %>% mutate(total=sum(n), proportion=n/total)
props$Migration<-as.factor(props$Migration)
ggplot(props, aes(x=zone_bin, y=proportion, fill=Migration))+
  geom_bar(position="stack", stat="identity")+
  facet_grid(category~season)





# are there are higher proportion of migrants in non-urban habitats?
mod.migration <- aov(n ~ zone_bin + season + category + Migration, data=richness_category)
summary(mod.migration)




################# Look at relationship between diet and habitat breadth and migratory status
habitat_migration <- merge(season_sp_habitat, avonet[,c(1,28)], by="SCIENTIFIC.NAME", by.y="Species2")
habitat_migration %>% group_by(Migration) %>% summarise(mean=mean(Habitat_breadth_IUCN))
# the widest habitat breadth is the partial migrants, which kinda makes sense
aov1 <- aov(log(Habitat_breadth_IUCN) ~ as.factor(Migration), data = habitat_migration)
ggplot(habitat_migration)+
  geom_boxplot(aes(x=as.factor(Migration), y=Habitat_breadth_IUCN))
TukeyHSD(aov1)


diet_migration <- merge(season_sp_diet, avonet[,c(1,28)], by="SCIENTIFIC.NAME", by.y="Species2")
diet_migration %>% group_by(Migration) %>% summarise(mean=mean(gini.index))
# full migrants have the most specialized diets
aov2 <- aov(gini.index ~ as.factor(Migration), data = diet_migration)

ggplot(diet_migration)+
  geom_boxplot(aes(x=as.factor(Migration), y=gini.index))
TukeyHSD(aov2)




#################### Looking at proportion of cells that are in urban vs not urban for each bird species
# see if that changes with season

season_uniquesp <- read.table("season_unique_species.txt", header=TRUE)
season_uniquesp <- season_uniquesp %>% 
  mutate(zone_bin = cut(abslat, breaks=c(0, 23.43621, 35, 50, 90), 
                        labels=c("Tropical", "Subtropical", "Temperate", "Subpolar")))
season_uniquesp <- season_uniquesp %>% 
  mutate(lat_bin = cut(abslat, breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80)))

# summarise: how many cells in each category for each season
summary.urban <- season_uniquesp %>% group_by(SCIENTIFIC.NAME, season, urban2) %>% count() %>% 
  group_by(SCIENTIFIC.NAME, season) %>% mutate(total=sum(n), prop=n/total)

summary.urban2 <- summary.urban %>% group_by(season, urban2) %>% summarise(mean=mean(prop))
# proportion in urban areas does not change with season, suggesting that birds are not moving into urban areas

# But I would only expect the proportion to increase at mid latitudes because that is where seasonality is highest
summary.urban <- season_uniquesp %>% group_by(SCIENTIFIC.NAME, season, urban2, zone_bin) %>% count() %>% 
  group_by(SCIENTIFIC.NAME, season, zone_bin) %>% mutate(total=sum(n), prop=n/total)

summary.urban2 <- summary.urban %>% group_by(season, urban2, zone_bin) %>% summarise(mean=mean(prop))

summary.urban2 %>% filter(urban2=="urban") %>% ggplot() +
  geom_point(aes(x=zone_bin, y=mean, color=season))+
  facet_wrap(~urban2)


# Look at finer resolution of latitude bins
summary.urban <- season_uniquesp %>% group_by(SCIENTIFIC.NAME, season, urban2, lat_bin) %>% count() %>% 
  group_by(SCIENTIFIC.NAME, season, lat_bin) %>% mutate(total=sum(n), prop=n/total)

summary.urban2 <- summary.urban %>% group_by(season, urban2, lat_bin) %>% summarise(mean=mean(prop))

summary.urban2 %>% filter(urban2=="urban") %>% ggplot() +
  geom_point(aes(x=lat_bin, y=mean, color=season))+
  geom_line(aes(x=lat_bin, y=mean, color=season, group=season))+
  facet_wrap(~urban2)







####### Now looking at mean latitude per species in summer and winter to see if natural birds are migrating more than urban
summary.lat <- season_uniquesp %>% group_by(SCIENTIFIC.NAME, season, urban2) %>% summarise(species.mean=mean(abslat))
# plot species level avg latitude as a boxplot
ggplot(summary.lat, aes(x=urban2, y=species.mean, fill=season))+
  geom_boxplot()


summary.lat2 <- summary.lat %>% group_by(season,urban2) %>% summarise(mean=mean(species.mean))

lat.anova <- aov(species.mean ~ season*urban2, dat=summary.lat)
summary(lat.anova)
library(marginaleffects)

means <- marginal_means(lat.anova, variables=c("season", "urban2"), cross=TRUE)

ggplot(means, aes(x=urban2, y=estimate, group=season, color=season))+
  geom_point(size=2, position=position_dodge(width=0.2))+
  geom_line(size=0.5, position=position_dodge(width=0.2))+
  scale_color_manual(values=c("grey30", "deepskyblue3"))+
  labs(y="Avg latitude")+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.15, position=position_dodge(width=0.2))









