############## Calculating threshold season ############
library(terra)
library(tidyverse)
library(sf)
library(iNEXT)
library(elevatr)

# here I will calculate the threshold I need seperately for summer and winter months using the iNEXT package

######## Summer ######

summer <- read.csv("summer_richness_summary.csv")
summer_top_cells <- summer %>% slice_max(total_SR, n=500) # find top 500 cells
unique(summer_top_cells$square)
summer_top_cells %>% group_by(square) %>% count()

names <- c("r1c1", "r1c2", "r1c3", "r1c4", "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4", "r3c2", "r3c3", "r3c4") # figure out names

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# filter out cells that I want
for (j in 4:4) {
  
  top_cell <- summer_top_cells %>% filter(square==names[j])
  # make list
  datalist = vector("list", length = length(years))
  
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[i], "_data/summer/", names[j], "_", years[i], "_summer_filt.txt", sep=""), header=TRUE, na.strings="")
    # filter out cells that are in the top cell
    dat_filt <- dat %>% filter(cell %in% top_cell$cell)
    dat_filt$cell <- as.integer(dat_filt$cell)
    # add to list
    datalist[[i]] <- dat_filt
  }
  # bind lists together
  data <- dplyr::bind_rows(datalist)
  # save as csv
  write.csv(data, paste("thresholding/summer/", names[j], "_topcells_sum.csv", sep=""))
  print(paste("finished", names[j]))
  rm(dat)
  rm(datalist)
  rm(data)
}


#########
# calculate sample size needed to reach certain coverage
coverage_list = vector("list", length = length(names))

for (j in 1:length(names)) {
  
  dat_top <- read.csv(paste("thresholding/summer/", names[j], "_topcells_sum.csv", sep=""))
  
  ebird.split <- dat_top %>% group_by(cell) %>% group_split()
  
  
  # make data frame for output
  output <- data.frame(matrix(NA,
                              nrow=length(ebird.split),
                              ncol=19))
  
  # make column names for output
  names(output)<-c("cell_ID", "SC_obs", "obs.richness", "sampsize_obs", "SC_80", "richness_80", "sampsize_80", "SC_90", "richness_90", 
                   "sampsize_90", "SC_95", 
                   "richness_95", "sampsize_95",
                   "SC_97", "richness_97", "sampsize_97", "SC_98", "richness_98", "sampsize_98")
  
  
  set.seed(20) # because it is bootstrapping
  withCallingHandlers ({
    for(i in 1:length(ebird.split)){
      output$cell_ID[i] <- mean(ebird.split[[i]]$cell)
      
      w <- length(warnings())
      # put into the format that is accepted by the iNEXT package
      temp <- ebird.split[[i]] %>% # select list element
        group_by(SCIENTIFIC.NAME, SAMPLING.EVENT.IDENTIFIER) %>% 
        summarize(present=n()) %>%
        mutate(present=1) %>%
        pivot_wider(names_from=SAMPLING.EVENT.IDENTIFIER, values_from=present, values_fill=0) %>%
        column_to_rownames(var="SCIENTIFIC.NAME") %>%
        as.data.frame() 
      
      # convert this dataframe into data format for iNext
      temp_inext <- as.incfreq(temp)
      #and now using estimateD to get qD
      out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=1000, nboot=50)
      out.inc.filt1<-out.inc$iNextEst$coverage_based %>% filter(Method=="Observed") # filter for observed coverage
      out.inc.filt2<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.80)==min(abs(SC-0.80)))  # filter for row that is closest to 95% coverage
      out.inc.filt3<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.90)==min(abs(SC-0.90)))  # filter for row that is closest to 95% coverage
      out.inc.filt4<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.95)==min(abs(SC-0.95)))  # filter for row that is closest to 95% coverage
      out.inc.filt5<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.97)==min(abs(SC-0.97))) # filter for row that is closest to 97% coverage
      out.inc.filt6<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.98)==min(abs(SC-0.98))) # filter for row that is closest to 98% coverage
      #sampling_profile <- estimateD(temp_inext, q=0, datatype="incidence_freq",base="coverage", 
      #                             level=0.80, conf=0.95) # set to a coverage level of 60
      
      # add output to a data frame
      output$SC_obs[i] <- out.inc.filt1$SC # sample coverage
      output$obs.richness[i] <- out.inc.filt1$qD
      output$sampsize_obs[i] <- out.inc.filt1$t
      output$SC_80[i] <- out.inc.filt2$SC # sample coverage
      output$richness_80[i] <- out.inc.filt2$qD
      output$sampsize_80[i] <- out.inc.filt2$t
      output$SC_90[i] <- out.inc.filt3$SC # sample coverage
      output$richness_90[i] <- out.inc.filt3$qD
      output$sampsize_90[i] <- out.inc.filt3$t
      output$SC_95[i] <- out.inc.filt4$SC # sample coverage
      output$richness_95[i] <- out.inc.filt4$qD
      output$sampsize_95[i] <- out.inc.filt4$t
      output$SC_97[i] <- out.inc.filt5$SC # sample coverage
      output$richness_97[i] <- out.inc.filt5$qD
      output$sampsize_97[i] <- out.inc.filt5$t
      output$SC_98[i] <- out.inc.filt6$SC # sample coverage
      output$richness_98[i] <- out.inc.filt6$qD
      output$sampsize_98[i] <- out.inc.filt6$t
      # add new columns to data frame
      print(paste("finished", i))
    }
  },
  warning = function(w){
    output$warnings[i] <<- w$message
    invokeRestart("muffleWarning")
  })
  
  output$square <- names[j]
  coverage_list[[j]] <- output
  print(paste("finished", names[j]))
}

coverage <- dplyr::bind_rows(coverage_list)

coverage %>% group_by(square) %>% summarise(n=n()) 

write.csv(coverage, "thresholding/summer/summer_coverage_top500.csv", row.names=FALSE)


################### WINTER #################
winter <- read.csv("winter_richness_summary.csv")
winter_top_cells <- winter %>% slice_max(total_SR, n=500) # find top 500 cells
unique(winter_top_cells$square)

names <- c("r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c3", "r2c4", "r3c2", "r3c3", "r3c4") # figure out names

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# filter out cells that I want
for (j in 9:length(names)) {
  
  top_cell <- winter_top_cells %>% filter(square==names[j])
  # make list
  datalist = vector("list", length = length(years))
  
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[i], "_data/winter/", names[j], "_", years[i], "_winter_filt.txt", sep=""), header=TRUE, na.strings="")
    # filter out cells that are in the top cell
    dat_filt <- dat %>% filter(cell %in% top_cell$cell)
    # add to list
    datalist[[i]] <- dat_filt
  }
  # bind lists together
  data <- dplyr::bind_rows(datalist)
  # save as csv
  write.csv(data, paste("thresholding/winter/", names[j], "_topcells_wint.csv", sep=""))
  print(paste("finished", names[j]))
  rm(dat)
  rm(datalist)
  rm(data)
}

#########
# calculate sample size needed to reach certain coverage
coverage_list = vector("list", length = length(names))

for (j in 5:length(names)) {
  
  dat_top <- read.csv(paste("thresholding/winter/", names[j], "_topcells_wint.csv", sep=""))
  
  ebird.split <- dat_top %>% group_by(cell) %>% group_split()
  
  
  # make data frame for output
  output <- data.frame(matrix(NA,
                              nrow=length(ebird.split),
                              ncol=19))
  
  # make column names for output
  names(output)<-c("cell_ID", "SC_obs", "obs.richness", "sampsize_obs", "SC_80", "richness_80", "sampsize_80", "SC_90", "richness_90", 
                   "sampsize_90", "SC_95", 
                   "richness_95", "sampsize_95",
                   "SC_97", "richness_97", "sampsize_97", "SC_98", "richness_98", "sampsize_98")
  
  
  set.seed(20) # because it is bootstrapping
  withCallingHandlers ({
    for(i in 1:length(ebird.split)){
      output$cell_ID[i] <- mean(ebird.split[[i]]$cell)
      
      w <- length(warnings())
      # put into the format that is accepted by the iNEXT package
      temp <- ebird.split[[i]] %>% # select list element
        group_by(SCIENTIFIC.NAME, SAMPLING.EVENT.IDENTIFIER) %>% 
        summarize(present=n()) %>%
        mutate(present=1) %>%
        pivot_wider(names_from=SAMPLING.EVENT.IDENTIFIER, values_from=present, values_fill=0) %>%
        column_to_rownames(var="SCIENTIFIC.NAME") %>%
        as.data.frame() 
      
      # convert this dataframe into data format for iNext
      temp_inext <- as.incfreq(temp)
      #and now using estimateD to get qD
      out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=1000, nboot=50)
      out.inc.filt1<-out.inc$iNextEst$coverage_based %>% filter(Method=="Observed") # filter for observed coverage
      out.inc.filt2<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.80)==min(abs(SC-0.80)))  # filter for row that is closest to 95% coverage
      out.inc.filt3<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.90)==min(abs(SC-0.90)))  # filter for row that is closest to 95% coverage
      out.inc.filt4<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.95)==min(abs(SC-0.95)))  # filter for row that is closest to 95% coverage
      out.inc.filt5<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.97)==min(abs(SC-0.97))) # filter for row that is closest to 97% coverage
      out.inc.filt6<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.98)==min(abs(SC-0.98))) # filter for row that is closest to 98% coverage
      #sampling_profile <- estimateD(temp_inext, q=0, datatype="incidence_freq",base="coverage", 
      #                             level=0.80, conf=0.95) # set to a coverage level of 60
      
      # add output to a data frame
      output$SC_obs[i] <- out.inc.filt1$SC # sample coverage
      output$obs.richness[i] <- out.inc.filt1$qD
      output$sampsize_obs[i] <- out.inc.filt1$t
      output$SC_80[i] <- out.inc.filt2$SC # sample coverage
      output$richness_80[i] <- out.inc.filt2$qD
      output$sampsize_80[i] <- out.inc.filt2$t
      output$SC_90[i] <- out.inc.filt3$SC # sample coverage
      output$richness_90[i] <- out.inc.filt3$qD
      output$sampsize_90[i] <- out.inc.filt3$t
      output$SC_95[i] <- out.inc.filt4$SC # sample coverage
      output$richness_95[i] <- out.inc.filt4$qD
      output$sampsize_95[i] <- out.inc.filt4$t
      output$SC_97[i] <- out.inc.filt5$SC # sample coverage
      output$richness_97[i] <- out.inc.filt5$qD
      output$sampsize_97[i] <- out.inc.filt5$t
      output$SC_98[i] <- out.inc.filt6$SC # sample coverage
      output$richness_98[i] <- out.inc.filt6$qD
      output$sampsize_98[i] <- out.inc.filt6$t
      # add new columns to data frame
      print(paste("finished", i))
    }
  },
  warning = function(w){
    output$warnings[i] <<- w$message
    invokeRestart("muffleWarning")
  })
  
  output$square <- names[j]
  coverage_list[[j]] <- output
  print(paste("finished", names[j]))
}

coverage <- dplyr::bind_rows(coverage_list)

coverage %>% group_by(square) %>% summarise(n=n()) 

write.csv(coverage, "thresholding/winter/winter_coverage_top500.csv", row.names=FALSE)


#############################
### Prepare data for modelling
sum_threshold <- read.csv("thresholding/summer/summer_coverage_top500.csv")

sum_dat <- read.csv("summer_richness_summary.csv")
hist(sum_threshold$obs.richness)
hist(sum_threshold$richness_80)
hist(sum_threshold$richness_90)
hist(sum_threshold$richness_95)
hist(sum_threshold$richness_97)
hist(sum_threshold$richness_98)

# look at 95th quantile
quantile(sum_threshold$sampsize_80, 0.95) # 21
quantile(sum_threshold$sampsize_90, 0.95) # 42
quantile(sum_threshold$sampsize_95, 0.95) # 86.5
quantile(sum_threshold$sampsize_97, 0.95) # 139.25
quantile(sum_threshold$sampsize_98, 0.95) # 201
# they are lower than the full year which is good

### Winter
wint_threshold <- read.csv("thresholding/winter/winter_coverage_top500.csv")

wint_dat <- read.csv("winter_richness_summary.csv")
hist(wint_threshold$obs.richness)
hist(wint_threshold$richness_80)
hist(wint_threshold$richness_90)
hist(wint_threshold$richness_95)
hist(wint_threshold$richness_97)
hist(wint_threshold$richness_98)

# look at 95th quantile
quantile(wint_threshold$sampsize_80, 0.95) # 34.5
quantile(wint_threshold$sampsize_90, 0.95) # 70
quantile(wint_threshold$sampsize_95, 0.95) # 120.5
quantile(wint_threshold$sampsize_97, 0.95) # 176.5
quantile(wint_threshold$sampsize_98, 0.95) # 249.5
# they are lower than the full year which is good, but not really that different. Interesting that winter is higher than summer.

# I will threshold them at their respective thresholds (might as well)
summer_filt95 <- sum_dat %>% filter(number_checklists >= 87) # 11109
winter_filt95 <- wint_dat %>% filter(number_checklists >= 87) # 28257
summary(wint_dat)

## Extract urbanization values for each
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtMollweide.tif")
# Summer
summer_filt95$x <- xFromCell(GHSL, summer_filt95$cell) # extract the coordinates from the cells
summer_filt95$y <- yFromCell(GHSL, summer_filt95$cell)

summer_filt95$urban <- as.data.frame(terra::extract(GHSL, summer_filt95[,c(9:10)]))$SMOD_global
summer_filt95$season <- "summer" # add season
summer_filt95 <- summer_filt95 %>% filter(!is.na(urban)) # remove NAs for urbanization - 20,611
# Winter
winter_filt95$x <- xFromCell(GHSL, winter_filt95$cell) # extract the coordinates from the cells
winter_filt95$y <- yFromCell(GHSL, winter_filt95$cell)

winter_filt95$urban <- as.data.frame(terra::extract(GHSL, winter_filt95[,c(9:10)]))$SMOD_global
winter_filt95$season <- "winter" # add season
winter_filt95 <- winter_filt95 %>% filter(!is.na(urban)) # remove NAs for urbanization - 21,658

## Put them together 
season_dat <- rbind(summer_filt95, winter_filt95)


sum_long <- summer_filt95 %>% select(cell, total_SR) %>% rename(summer=total_SR)
wint_long <- winter_filt95 %>% select(cell, total_SR) %>% rename(winter=total_SR)

season_long <- merge(sum_long, wint_long, by="cell")
season_long$x <- xFromCell(GHSL, season_long$cell)
season_long$y <- yFromCell(GHSL, season_long$cell)

season_long <- season_long %>% mutate(diff = summer - winter)
world <- ne_countries(scale = "medium", returnclass = "sf")
# plot
full.map<-ggplot(data=world)+
  geom_sf(lwd=0.15, fill="white") +
  geom_point(data=season_long, aes(x=x, y=y, color=diff, size=diff), size=.65, alpha=0.5) +
  #coord_sf(crs = 4326, expand = TRUE, xlim=c(-180, 180), ylim=c(-55, 70)) +
  labs(x="Longitude", y="Latitude")+
  theme_classic()+
  scale_color_viridis_c()+
  # geom_hline(yintercept=c(23.4, -23.4, 35, -35, 50, -50, 66.5, -66.5), alpha=0.7, lty=3)+ # geographic zones
  geom_hline(yintercept=0, alpha=0.8, lty=3)  # for equator
#  theme_classic() + theme(legend.title=element_blank(), legend.position = "none", text=element_text(size=15),
 #                         axis.text = element_blank(), axis.ticks = element_blank(), 
  #                        axis.title=element_blank(), axis.line = element_blank(), strip.background=element_blank(), strip.text=element_blank()) 
full.map
## Assign hemisphere
#season_dat <- season_dat %>% mutate(hemisphere = if_else(x>0, "northern", "southern"))

## Add continent
continents <- st_read("/Volumes/Backup/eBird/continent-poly/Continents.shp")
continents_moll <- st_transform(continents, crs=crs(GHSL)) 

season_dat_sf <- st_as_sf(season_dat, coords=c("x", "y"), crs=st_crs(GHSL))

dat_cont <- st_join(season_dat_sf, continents_moll[,"CONTINENT"], left=TRUE, join=st_nearest_feature) # joining by nearest feature

## Extract biome
biomes <- st_read("/Volumes/Backup/eBird/wwf_biomes/wwf_terr_ecos.shp")
biomes_moll <- st_transform(biomes, crs=crs(GHSL)) 
dat_withbiome <- st_join(dat_cont, biomes_moll[,c("REALM", "BIOME")], left=TRUE, join=st_nearest_feature)

## Extract elevation
dat_latlong <- st_transform(dat_withbiome, crs=st_crs(4326)) # get lat long coordinates as well for the elevation extraction
latlong_df <- dat_latlong %>% mutate(long = sf::st_coordinates(.)[,1],
                                     lat = sf::st_coordinates(.)[,2])

latlong_df <- get_elev_point(latlong_df[,c(15,16,2:13)], prj=crs(dat_latlong), src="aws", overwrite=TRUE) # extract elevations from amazon web services

dat <- as.data.frame(st_transform(latlong_df, crs=crs(GHSL)) %>% mutate(x = sf::st_coordinates(.)[,1],
                                                                             y = sf::st_coordinates(.)[,2]))


######### Continue to prepare data for modelling

# assign hemisphere
dat$hemisphere <- "northern"
dat$hemisphere[dat$lat<0]<-"southern"

dat$urban<-as.factor(dat$urban) # make urban score a factor (instead of numeric)
dat$BIOME <- as.factor(dat$BIOME) # make biome a factor

# make another column with just 3 categories
dat <- dat %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat$urban2 <- as.factor(dat$urban2)
dat$abslat <- abs(dat$lat) # absolute latitude

# Divide by quartiles
#for (i in 1:nrow(dat)){
#  if (dat$long[i] < 0 & dat$hemisphere[i] == "northern") { # quadrant 1 is North America
#    dat$quadrant[i] <- 1
#  }
#  else if (dat$long[i] > 0 & dat$hemisphere[i] == "northern") { # quadrant 2 is europe and asia and N Africa
#    dat$quadrant[i] <- 2
#  }
#  else if (dat$long[i] < 0 & dat$hemisphere[i] == "southern") { # quadrant 3 is south america 
#    dat$quadrant[i] <- 3 
#  }
#  else {dat$quadrant[i] <- 4} # quadrant 4 is oceania and southern africa
#}
#dat$quadrant <- as.factor(dat$quadrant)

#dat <- dat %>% filter(!CONTINENT == "Antarctica") # filter out antarctica


urb <- dat %>% filter(urban2==3)
range(urb$lat) # the highest latitude is 64.15 and the lowest is -54.84

suburb <- dat %>% filter(urban2==2)
range(suburb$lat) # highest latitude is 65.76 and lowest latitude is -53.11

nat <- dat %>% filter(urban2==1)
range(nat$lat) # highest latitude is 65.76 and lowest latitude is -53.11
# 42266 - 42243
dat <- dat %>% filter(lat <= 66 & lat >= -55) # cut off data at the highest latitude range

dat$season <- factor(dat$season, levels = c("winter", "summer"),
                     labels = c("Winter", "Summer"))

dat$urban2 <- factor(dat$urban2, levels = c("1", "2", "3"),
                     labels = c("Natural", "Suburban", "Urban"))


#### Add precipitation data
precip <- rast("precipitation/wc2.1_5m_bio_12.tif")
#dat <- read.csv("season_modeling_data.csv")
dat$precip <- as.data.frame(terra::extract(precip, dat[,c(1:2)], method="bilinear"))$wc2.1_5m_bio_12

## Save data
write_csv(dat, "season_modeling_data.csv")








