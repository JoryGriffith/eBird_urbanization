############### Code to scale up the analysis to 5x5 km rasters as a sensitivity analysis #####################

library(terra)
library(tidyverse)
library(lubridate)
library(beepr)
library(sf)
library(foreach)
library(doParallel)
library(iNEXT)


# load original reprojected raster file

#GHSL <- rast("/Volumes/Expansion/eBird/SMOD_global/reprojected.SMOD_global.tif")

#GHSL_5km <- aggregate(GHSL, fact=5, fun="mean", cores=4, filename="/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_cellsize.tif")

#plot(GHSL_5km)

GHSL_5km <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_cellsize.tif")
# yay! created the new raster

# Now I need to aggregate the eBird data into that raster and summarise

########## 1) extract cell no. for each point ########
names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# trying to run in parallel
n.cores <- parallel::detectCores()

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()


for (j in 6:6) {

# foreach  (i = 6:length(names),
 #           .combine = 'c') %dopar% {
for(i in 6:length(names)){
  dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_", years[j],"_data/custom_bbox/", names[i], "_", years[j], "_filt.txt", sep=""), 
                    header=TRUE, na.strings="")
  # turn into spatvector
  
  vect <- vect(dat, crs=crs(GHSL_5km),geom=c("LONGITUDE","LATITUDE"))
  
  xy=geom(vect)
  
  # get cell number that each point is in
  dat$cell_5km_scale<-cellFromXY(GHSL_5km, xy[,3:4])
  # also get coordinates for the midpoint of each cell
  write.table(dat, paste("/Volumes/Expansion/eBird/eBird_", years[j],"_data/custom_bbox/", names[i], "_", years[j], "_filt.txt", sep=""), row.names=FALSE)
  print(paste("finished", names[i]))
  rm(dat)
}
print(paste("finished", years[j]))
}
beep()
beep_on_error()




############# 2) Summarise by cell no. ########
for (j in 5:9){ # skipped 5-9, need to run on other computer
  datalist = vector("list", length = length(years))
  # loop for each year
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_", years[i], "_data/custom_bbox/", names[j], "_", years[i], "_filt.txt", sep=""), 
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
  
  # summarise
  dat2$month <- month(dat2$OBSERVATION.DATE)
  
  # make dataframe of months and number of checklists
  months <- dat2 %>% group_by(cell_5km_scale, month) %>% 
    summarise(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)))
  
  months_wide <- pivot_wider(months, id_cols = cell_5km_scale, names_from=month, values_from = number_checklists)
  months_wide[is.na(months_wide)]<-0
  
  ## aggregate to get number of checklists and raw richness per cell
  dat_summary <- dat2 %>%
    group_by(cell_5km_scale) %>%
    summarize(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)),
              total_SR=length(unique(SCIENTIFIC.NAME)),
              total_duration=sum(DURATION.MINUTES),
              avg_duration=mean(DURATION.MINUTES),
              no_months=length(unique(month))
    )
  
  dat_summary2 <- merge(dat_summary, months_wide, by="cell_5km_scale")
  dat_summary2$square=names[j]
  # save as csv
  write.csv(dat_summary2, paste("5yr_summary_5km/", names[j], "_SR_5km.csv", sep=""))
  print(paste("finished", names[j]))
}

##### put all together
list_csv_files <- list.files(path = "5yr_summary_5km/", pattern="*.csv")

names <- tolower(gsub('_SR_5km.csv', "", list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary_5km/", list_csv_files[i], sep="")))
}

dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)
#save all the summaries as a csv
write.csv(dat, "global_richness_summary_5km.csv", row.names=FALSE)

# find top 500 cells
top_cells <- dat %>% slice_max(total_SR, n=500)
write.csv(top_cells, "top_500_cells_5km.csv", row.names=FALSE)


############# 3) Calculate threshold for inclusion ##############
top_cells_5km <- read.csv("top_500_cells_5km.csv") 

top_cells_5km %>% group_by(square) %>% summarise(n=n())
unique(top_cells_5km$square)

names <- c("r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c3", "r2c4", "r3c2", "r3c3")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# make this in parallel

# filter out the cells that I want
foreach  (j = 1:length(names),
                     .combine = 'c') %dopar% {

#for (j in 1:length(names)) {
  
  top_cell <- top_cells_5km %>% filter(square==names[j])
  # make list
  datalist = vector("list", length = length(years))
  
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_", years[i], "_data/custom_bbox/", names[j], "_", years[i], "_filt.txt", sep=""), header=TRUE)
    # filter out cells that are in the top cell
    dat_filt <- dat %>% filter(cell_5km_scale %in% top_cell$cell_5km_scale)
    # add to list
    datalist[[i]] <- dat_filt
  }
  # bind lists together
  data <- dplyr::bind_rows(datalist)
  # save as csv
  write.csv(data, paste("thresholding_5km/", names[j], "_topcells_5km.csv", sep=""), row.names=FALSE)
  print(paste("finished", names[j]))
  rm(dat)
  rm(datalist)
  rm(data)
}


### Make species accumulation curves and calculate coverage
coverage_list = vector("list", length = length(names))


# do in parallel
#foreach  (j = 1:length(names),
 #           .combine = 'c') %dopar% {
for (j in 2:length(names)){ # 73 in r2c1
  
  dat_top <- read.csv(paste("thresholding_5km/", names[j], "_topcells_5km.csv", sep=""))
  
  ebird.split <- dat_top %>% group_by(cell_5km_scale) %>% group_split()
  
  
  # make data frame for output
  output <- data.frame(matrix(NA,
                              nrow=length(ebird.split),
                              ncol=13))
  
  # make column names for output
  names(output)<-c("cell_ID", "SC_obs", "obs.richness", "sampsize_obs", "SC_95", "richness_95", "sampsize_95",
                   "SC_97", "richness_97", "sampsize_97", "SC_98", "richness_98", "sampsize_98")
  
  
  set.seed(20) # because it is bootstrapping
  withCallingHandlers ({
    for(i in 1:length(ebird.split)){
      output$cell_ID[i] <- mean(ebird.split[[i]]$cell_5km_scale)
      
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
      out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=500, nboot=50)
      out.inc.filt1<-out.inc$iNextEst$coverage_based %>% filter(Method=="Observed") # filter for observed coverage
      out.inc.filt2<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.95)==min(abs(SC-0.95)))  # filter for row that is closest to 95% coverage
      out.inc.filt3<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.97)==min(abs(SC-0.97))) # filter for row that is closest to 97% coverage
      out.inc.filt4<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.98)==min(abs(SC-0.98))) # filter for row that is closest to 98% coverage
      #sampling_profile <- estimateD(temp_inext, q=0, datatype="incidence_freq",base="coverage", 
      #                             level=0.80, conf=0.95) # set to a coverage level of 60
      
      # add output to a data frame
      output$SC_obs[i] <- out.inc.filt1$SC # sample coverage
      output$obs.richness[i] <- out.inc.filt1$qD
      output$sampsize_obs[i] <- out.inc.filt1$t
      output$SC_95[i] <- out.inc.filt2$SC # sample coverage
      output$richness_95[i] <- out.inc.filt2$qD
      output$sampsize_95[i] <- out.inc.filt2$t
      output$SC_97[i] <- out.inc.filt3$SC # sample coverage
      output$richness_97[i] <- out.inc.filt3$qD
      output$sampsize_97[i] <- out.inc.filt3$t
      output$SC_98[i] <- out.inc.filt4$SC # sample coverage
      output$richness_98[i] <- out.inc.filt4$qD
      output$sampsize_98[i] <- out.inc.filt4$t
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

write.csv(coverage, "coverage_top500_5km.csv", row.names=FALSE)


#################################
#################################

###### Scale up raster for season

#### Extract cell no. for season
GHSL_5km <- rast("/Volumes/Expansion/eBird/SMOD_global/SMOD_5km_cellsize.tif")

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

for (j in 1:length(years)){
  for (i in 1:length(names)){
    
    dat <- read.delim(paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), header=TRUE, na.strings="")
    # turn into spatvector
    vect <- vect(dat, crs=crs(GHSL_5km),geom=c("LONGITUDE","LATITUDE"))
    
    xy=geom(vect)
    
    # get cell number that each point is in
    dat$cell_5km<-cellFromXY(GHSL_5km, xy[,3:4])
    # also get coordinates for the midpoint of each cell
    write.table(dat, paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
    
  }
  print(paste("finished", years[j]))
}


# Winter
for (j in 1:length(years)){
  for (i in 1:length(names)){
    
    dat <- read.delim(paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), header=TRUE, na.strings="")
    # turn into spatvector
    vect <- vect(dat, crs=crs(GHSL_5km),geom=c("LONGITUDE","LATITUDE"))
    
    xy=geom(vect)
    
    # get cell number that each point is in
    dat$cell_5km<-cellFromXY(GHSL_5km, xy[,3:4])
    # also get coordinates for the midpoint of each cell
    write.table(dat, paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""), row.names=FALSE)
    print(paste("finished", names[i]))
    
  }
  print(paste("finished", years[j]))
}


#### Summarize by season ########
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
    dat <- read.table(paste("/Volumes/Expansion/eBird/eBird_", years[i], "_data/summer/", names[j], "_", years[i], "_summer_filt.txt", sep=""), 
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
  write.csv(dat_summary, paste("5yr_summary_5km/summer/", names[j], "_summer_SR_5km.csv", sep=""))
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
  write.csv(dat_summary, paste("5yr_summary_5km/winter/", names[j], "_winter_SR_5km.csv", sep=""))
  print(paste("finished", names[j]))
}

# Put all together - winter
list_csv_files <- list.files(path = "5yr_summary_5km/winter/", pattern="*.csv")
#dat <- readr::read_csv(paste("5yr_summary/", list_csv_files, sep=""), id = "file_name")

names <- tolower(gsub('_winter_SR_5km.csv', "", list_csv_files))

for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary_5km/winter/", list_csv_files[i], sep="")))
}
class(r4c3$square)
r4c1$square <- as.character(r4c1$square)
r4c3$square <- as.character(r4c3$square)
dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)
write.csv(dat, "winter_richness_summary_5km.csv", row.names=FALSE)

# Put all together - summer
list_csv_files <- list.files(path = "5yr_summary_5km/summer/", pattern="*.csv")
#dat <- readr::read_csv(paste("5yr_summary/", list_csv_files, sep=""), id = "file_name")

names <- tolower(gsub('_summer_SR_5km.csv', "", list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary_5km/summer/", list_csv_files[i], sep="")))
}

dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)

#save all the summaries as a csv
write.csv(dat, "summer_richness_summary_5km.csv", row.names=FALSE)

























