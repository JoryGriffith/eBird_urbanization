################## Summarise each cell in each season ########

library(tidyverse)
library(lubridate)
library(terra)
library(sf)

######### 
years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

## Summer
for (j in 1:length(names)){
  datalist = vector("list", length = length(years))
  # loop for each year
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[i], "_data/summer/", names[j], "_", years[i], "_summer_filt.txt", sep=""), 
                      header=TRUE) # load data
    
    dat$SCIENTIFIC.NAME <- as.character(dat$SCIENTIFIC.NAME)
    dat$OBSERVATION.DATE <- as.character(dat$OBSERVATION.DATE)
    dat$OBSERVER.ID <- as.character(dat$OBSERVER.ID)
    dat$SAMPLING.EVENT.IDENTIFIER <- as.character(dat$SAMPLING.EVENT.IDENTIFIER)
    dat$OBSERVATION.COUNT <- as.character(dat$OBSERVATION.COUNT)
    dat$GROUP.IDENTIFIER <- as.character(dat$GROUP.IDENTIFIER)
    datalist[[i]] <- dat # put in a list
  }
  # bind lists together
  
  dat2 <- dplyr::bind_rows(datalist)
  dat2$month <- month(dat2$OBSERVATION.DATE)
  # summarise
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat2 %>%
    group_by(cell) %>%
    summarize(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)),
              total_SR=length(unique(SCIENTIFIC.NAME)),
              total_duration=sum(DURATION.MINUTES),
              avg_duration=mean(DURATION.MINUTES),
              no_months=length(unique(month))
    )
  dat_summary$square=names[j]
  # save as csv
  write.csv(dat_summary, paste("5yr_summary/summer/", names[j], "_summer_SR.csv", sep=""))
  print(paste("finished", names[j]))
}


## Winter
for (j in 1:length(names)){
  datalist = vector("list", length = length(years))
  # loop for each year
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[i], "_data/winter/", names[j], "_", years[i], "_winter_filt.txt", sep=""), 
                      header=TRUE) # load data
    
    dat$SCIENTIFIC.NAME <- as.character(dat$SCIENTIFIC.NAME)
    dat$OBSERVATION.DATE <- as.character(dat$OBSERVATION.DATE)
    dat$OBSERVER.ID <- as.character(dat$OBSERVER.ID)
    dat$SAMPLING.EVENT.IDENTIFIER <- as.character(dat$SAMPLING.EVENT.IDENTIFIER)
    dat$OBSERVATION.COUNT <- as.character(dat$OBSERVATION.COUNT)
    dat$GROUP.IDENTIFIER <- as.character(dat$GROUP.IDENTIFIER)
    datalist[[i]] <- dat # put in a list
  }
  # bind lists together
  
  dat2 <- dplyr::bind_rows(datalist)
  dat2$month <- month(dat2$OBSERVATION.DATE)
  # summarise
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat2 %>%
    group_by(cell) %>%
    summarize(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)),
              total_SR=length(unique(SCIENTIFIC.NAME)),
              total_duration=sum(DURATION.MINUTES),
              avg_duration=mean(DURATION.MINUTES),
              no_months=length(unique(month))
    )
  dat_summary$square=names[j]
  # save as csv
  write.csv(dat_summary, paste("5yr_summary/winter/", names[j], "_winter_SR.csv", sep=""))
  print(paste("finished", names[j]))
}

############ put all together - winter
list_csv_files <- list.files(path = "5yr_summary/winter/", pattern="*.csv")
#dat <- readr::read_csv(paste("5yr_summary/", list_csv_files, sep=""), id = "file_name")

names <- tolower(gsub('_winter_SR.csv', "", list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary/winter/", list_csv_files[i], sep="")))
}
class(r4c3$square)
r4c1$square <- as.character(r4c1$square)
r4c3$square <- as.character(r4c3$square)
dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)


#save all the summaries as a csv
write.csv(dat, "winter_richness_summary.csv", row.names=FALSE)


############ put all together - summer
list_csv_files <- list.files(path = "5yr_summary/summer/", pattern="*.csv")
#dat <- readr::read_csv(paste("5yr_summary/", list_csv_files, sep=""), id = "file_name")

names <- tolower(gsub('_summer_SR.csv', "", list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary/summer/", list_csv_files[i], sep="")))
}

dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)

#save all the summaries as a csv
write.csv(dat, "summer_richness_summary.csv", row.names=FALSE)


#### Add the urbanization data
summer <- read.csv("summer_richness_summary.csv")
winter <- read.csv("winter_richness_summary.csv")
# similar number of data points

# extract urbanization values for each
# summer
GHSL <- rast("/Volumes/Backup/eBird/SMOD_global/GHSL_filtered.tif")

summer$x <- xFromCell(GHSL, summer$cell) # extract the coordinates from the cells
summer$y <- yFromCell(GHSL, summer$cell)

summer$urban <- as.data.frame(terra::extract(GHSL, summer[,c(9:10)]))$SMOD_global
summer$season <- "summer" # add season

# winter
winter$x <- xFromCell(GHSL, winter$cell) # extract the coordinates from the cells
winter$y <- yFromCell(GHSL, winter$cell)

winter$urban <- as.data.frame(terra::extract(GHSL, winter[,c(9:10)]))$SMOD_global
winter$season <- "winter"

# plot coverage of both
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")

summer_coverage_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=summer, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude", color="SR")+
  theme_bw()

winter_coverage_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=winter, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude", color="SR")+
  theme_bw()

# put them together
season_dat <- rbind(summer, winter)
# almost 2 million points

# threshold to 87 (like the other models)
season_dat_filt <- season_dat %>% filter(number_checklists >=87)
season_dat_filt <- season_dat_filt %>% filter(!is.na(urban)) # remove NAs
# brings it down to 37,058
season_dat_filt %>% group_by(season) %>% summarise(n=n()) # 16,116 summer and 20,942 winter

# plot by season
seasonal_coverage_plot <- ggplot(data=world)+
  geom_sf() +
  geom_point(data=season_dat_filt, aes(x=x, y=y, color=total_SR), size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude", color="SR")+
  theme_bw() +
  facet_wrap(~season)
# pretty similar coverage, that's good

season_dat_filt %>% group_by(urban) %>% summarise(n=n())

urb_levels <- season_dat_filt

urb_levels$urban[which(urb_levels$urban %in% c(11, 12, 13))] <- 1
urb_levels$urban[which(urb_levels$urban %in% c(21, 22, 23))] <- 2
urb_levels$urban[which(urb_levels$urban==30)] <- 3

urb_levels$urban <- as.factor(urb_levels$urban)
urban_names <- c(
  `1` = "Rural",
  `2` = "Peri-Urban",
  `3` = "Urban"
)

LDG <- ggplot(urb_levels, aes(x=abs(y), y=total_SR, color=urban)) +
  geom_point(alpha=0.1, shape=1)+
  geom_smooth(method="lm")+
  labs(x="Absolute Latitude", y="Species Richness")+
  theme_bw() +
  facet_wrap(~season)# has the expected relationship!

##################
# extract continent 
continents <- st_read("/Volumes/Backup/eBird/continent-poly/Continents.shp")
#plot(continents)

dat_sf <- st_as_sf(season_dat_filt, coords=c('x', "y"), crs=st_crs(continents))

dat_cont <- st_join(dat_sf, continents[,"CONTINENT"], left=TRUE, join=st_nearest_feature) # joining by nearest feature

###################
# extract biome 
biomes <- st_read("/Volumes/Backup/eBird/wwf_biomes/wwf_terr_ecos.shp")

dat_withbiome <- st_join(dat_cont, biomes[,"BIOME"], left=TRUE, join=st_nearest_feature)

# create seperate columns for lat long again
datFINAL <- as.data.frame(dat_withbiome[,-1] %>% mutate(long = sf::st_coordinates(.)[,1],
                                                        lat = sf::st_coordinates(.)[,2]))

datFINAL <- datFINAL %>% mutate(hemisphere = if_else(lat>0, "northern", "southern"))

# save final data as csv
write_csv(datFINAL, "season_model_data.csv")








