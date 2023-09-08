###### This is where I use iNEXT to figure out the threshold that I need to reach 95% in the species accumulation curves

library(terra)
library(tidyverse)
library(sf)
library(beepr)
library(scales)
library(rnaturalearth)
library(rnaturalearthdata)

# First I want to load and look at the top cells so that I know where they are coming from
top_cells <- read.csv("top_500_cells.csv") 
# there are 504 because some of them are tied
top_cells %>% group_by(square) %>% summarise(n=n())
unique(top_cells$square) # figure out which squares the top cells are in

#########
# Need to extract the raw data from these cells and put it into a data frame so I can calculate coverage

names <- c("r1c1", "r1c2", "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4", "r3c2", "r3c3")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# make loop for filter out cells that I want
for (j in 1:length(names)) {
  
  top_cell <- top_cells %>% filter(square==names[j])
# make list
   datalist = vector("list", length = length(years))

for (i in 1:length(years)) {
  dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_", years[i], "_data/custom_bbox/", names[j], "_", years[i], "_filt.txt", sep=""), header=TRUE)
  # filter out cells that are in the top cell
  dat_filt <- dat %>% filter(cell %in% top_cell$cell)
  # add to list
  datalist[[i]] <- dat_filt
  }
# bind lists together
  data <- dplyr::bind_rows(datalist)
# save as csv
write.csv(data, paste("thresholding/", names[j], "_topcells.csv", sep=""))
print(paste("finished", names[j]))
rm(dat)
rm(datalist)
rm(data)
}



#######################
# Calculate coverage and sample size needed to reach 95%

# make datalist for output
coverage_list = vector("list", length = length(names))

for (j in 1:length(names)) {
  
  dat_top <- read.csv(paste("thresholding/", names[j], "_topcells.csv", sep=""))
  
  ebird.split <- dat_top %>% group_by(cell) %>% group_split()
  
  
  # make data frame for output
  output <- data.frame(matrix(NA,
                              nrow=length(ebird.split),
                              ncol=13))
  
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

write.csv(coverage, "thresholding/coverage_top500.csv", row.names=FALSE)



###############################
# Look at coverage threshold
coverage <- read.csv("thresholding/coverage_top500.csv")
summary(coverage)
# look at distribution of sample sizes
hist(coverage$obs.richness)
hist(coverage$richness_95)
hist(coverage$richness_97)
hist(coverage$richness_98)

# look at mean richness at different coverages
mean(coverage$obs.richness) # 363.9464
mean(coverage$richness_95) # 204.684
mean(coverage$richness_97) # 232.3299
mean(coverage$richness_98) # 254

# look at the mean of sample size at different coverages
mean(coverage$sampsize_obs) # 2310.851
mean(coverage$sampsize_95) # 46
mean(coverage$sampsize_97) # 75
mean(coverage$sampsize_98) # 111

# look at 95th quantile
quantile(coverage$sampsize_95, 0.95) # 87
quantile(coverage$sampsize_97, 0.95) # 139
quantile(coverage$sampsize_98, 0.95) # 203



##############################
# Apply different thresholds to the data

# load summary data
summary <- read.csv("global_richness_summary.csv") # 2.1 mil

# threshold for 95 coverage
summary_filt95 <- summary %>% filter(number_checklists >= 87) # 79,768
# threshold for 97 coverage
summary_filt97 <- summary %>% filter(number_checklists >= 139) # 53,066
# threshold for 98 coverage
summary_filt98 <- summary %>% filter(number_checklists >= 203) # 37,081

# Going to use the 95% coverage filter for now

######## Add in lat long and urbanization score


# now need to put this in raster form and merge with the urbanization raster
# this raster layer has the places with high human impact removed
GHSL<-rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif.tif")

# need to extract cell numbers, urbanization scores, and x and y coordinates from the raster
summary_filt95$x <- xFromCell(GHSL, summary_filt95$cell) # extract the coordinates from the cells
summary_filt95$y <- yFromCell(GHSL, summary_filt95$cell)

summary_filt95$urban <- as.data.frame(terra::extract(GHSL, summary_filt95[,c(21:22)]))$SMOD_global

summary_filt95 %>% group_by(urban) %>% summarise(n=n())

# remove ones with NaN urbanization score
summary_filt <- summary_filt95 %>% filter(!is.nan(urban)) # 70749 datapoints now
# save the thresholded data
write.csv(summary_filt, "5yr_summary/summary_thresholded.csv", row.names=FALSE)


###############################################################################################
#### Prepare data for modelling

world <- ne_countries(scale = "medium", returnclass = "sf")

#world <- ne_countries(scale = "medium", type="map_units", returnclass = "sf")
GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtered.tif")

# load thresholded summary data
dat <- read.csv("5yr_summary/summary_thresholded.csv")
dat <- rename(dat, lat = y, long = x)

# First I need to get the other data that I want to include in the model

######## Assign hemisphere
dat <- dat %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))

######### Assign continent
#dat_sf <- st_as_sf(dat, coords=c('long', "lat"), crs=st_crs(world))

#joined <- st_join(dat_sf, world)


#dat_joined <- as.data.frame(joined[,c(1:22,41)] %>% mutate(long = sf::st_coordinates(.)[,1],
#         lat = sf::st_coordinates(.)[,2]))# just keep country name
# extract continent using country name
#dat_joined$continent <- countrycode(sourcevar = dat_joined[,"name_long"],
#                                    origin = "country.name",
#                                     destination = "continent")


##################
# Trying new way to do continent
continents <- st_read("/Volumes/Expansion/eBird/continent-poly/Continents.shp")
#plot(continents)

dat_sf <- st_as_sf(dat, coords=c('long', "lat"), crs=st_crs(continents))

dat_cont <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature) # joining by nearest feature


#########################
# Extract biome
# Classifying points into biomes using a terrestrial biomes shapefile
biomes <- st_read("/Volumes/Expansion/eBird/wwf_biomes/wwf_terr_ecos.shp")
class(biomes) # sf and data frame
# look at how many biomes there are
length(unique(biomes$REALM)) # 9 of these
length(unique(biomes$BIOME)) # 16 biomes
length(unique(biomes$ECO_NAME)) # 827 ecoregion names

# plot biome
#plot(biomes["BIOME"])

# want to extract biomes
dat_withbiome <- st_join(dat_cont, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)


# create seperate columns for lat long again
datFINAL <- as.data.frame(dat_withbiome[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
                                                        lat = sf::st_coordinates(.)[,2]))

summary(datFINAL)
# save as csv

#write_csv(datFINAL, "modeling_data.csv")

################################
## Add elevation
#dat.mod <- read_csv("modeling_data.csv")
dat.mod2 <- datFINAL[,c(25,26,1:24)] # reorder because lat and long need to be the first and second column

dat.mod.ele <- get_elev_point(dat.mod2, prj=crs(GHSL), src="aws") # extract elevations from amazon web services

dat.mod.ele.df <- as.data.frame(dat.mod.ele) %>% rename(long=coords.x1, lat=coords.x2)

write_csv(dat, "modeling_data.csv")








