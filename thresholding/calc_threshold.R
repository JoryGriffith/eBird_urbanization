###### This is where I use iNEXT to figure out the threshold that I need to reach 95% in the species accumulation curves

library(terra)
library(tidyverse)
library(sf)
library(auk)
library(iNEXT)
library(beepr)
library(scales)

# First I want to load and look at the top cells so that I know where they are coming from
top_cells <- read.csv("5yr_summary/top_500_cells.csv") # this is the updated top cells that count across years
# there are 504 because some of them are tied

unique(top_cells$square)
# "r3c2"    "r2c2ABA" "r2c2AA"  "r2c2ABB" "r2c3"    "r2c1"    "r2c4"    "r2c2B"   
# they are all in the 3rd or 2ns row which is good, means they are close to the equator

top_cells %>% group_by(square) %>% summarise(n=n())
# figure out which squares the top cells are in

#########
# Need to extract the raw data from these cells and put it into a data frame so I can calculate coverage

names <- c("r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4", "r3c2")
years <- c(2017, 2018, 2019, 2020, 2021)




##### r2c1
for (j in 3:length(names)) {
  
  top_cell <- top_cells %>% filter(square==names[j])
# make list
  datalist = vector("list", length = length(years))

for (i in 1:length(years)) {
  dat <- read.table(paste("eBird_", years[i], "_data/custom_bbox/", names[j], "_", years[i], "_filt.txt", sep=""), header=TRUE)
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
}


#######################
# Calculate coverage and sample size needed to reach 95%

# make datalist for output
coverage_list = vector("list", length = length(names))

for (j in 1:length(names)) {
  
  dat_top <- read.csv(paste("thresholding/", names[j], "_topcells.csv", sep=""))
  
  ebird.split <- dat_top %>% group_by(cell) %>% group_split()
  # this is a bit unnecessary because this is only one cell
  
  # make data frame for output
  output <- data.frame(matrix(NA,
                              nrow=length(ebird.split),
                              ncol=7))
  
  # make column names for output
  names(output)<-c("cell_ID", "SC_obs", "obs.richness", "sampsize_obs", "SC_95", "richness_95", "sampsize_95")
  
  
  set.seed(20) # because it is bootstrapping
  withCallingHandlers ({
    for(i in 1:length(ebird.split)){
      output$cell_ID[i] <- mean(ebird.split[[i]]$cell)
      
      w <- length(warnings())
      # put into the format that is accepted by the iNEXT package
      temp <- ebird.split[[i]] %>% # select list element
        group_by(scientific_name, sampling_event_identifier) %>% 
        summarize(present=n()) %>%
        mutate(present=1) %>%
        pivot_wider(names_from=sampling_event_identifier, values_from=present, values_fill=0) %>%
        column_to_rownames(var="scientific_name") %>%
        as.data.frame() 
      
      # convert this dataframe into data format for iNext
      temp_inext <- as.incfreq(temp)
      #and now using estimateD to get qD
      out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=500, nboot=50)
      out.inc.filt1<-out.inc$iNextEst$coverage_based %>% filter(Method=="Observed") # filter for observed coverage
      out.inc.filt2<-out.inc$iNextEst$coverage_based %>% filter(abs(SC-0.95)==min(abs(SC-0.95)))  # filter for row that is closest to 95% coverage
      
      #sampling_profile <- estimateD(temp_inext, q=0, datatype="incidence_freq",base="coverage", 
      #                             level=0.80, conf=0.95) # set to a coverage level of 60
      
      # add output to a data frame
      output$SC_obs[i] <- out.inc.filt1$SC # sample coverage
      output$obs.richness[i] <- out.inc.filt1$qD
      output$sampsize_obs[i] <- out.inc.filt1$t
      output$SC_95[i] <- out.inc.filt2$SC # sample coverage
      output$richness_95[i] <- out.inc.filt2$qD
      output$sampsize_95[i] <- out.inc.filt2$t
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

write.csv(coverage, "thresholding/coverage_top500.csv", row.names=FALSE)


###############################
# Look at coverage threshold
coverage <- read.csv("thresholding/coverage_top500.csv")

# look at distribution of sample sizes
hist(coverage$obs.richness)
mean(coverage$obs.richness)
hist(coverage$richness_95)
mean(coverage$richness_95) # definitely very different from the observed richness

mean(coverage$sampsize_obs) # 752 observed sample size
mean(coverage$sampsize_95) # 53 mean sample size for 95% coverage
quantile(coverage$sampsize_95, 0.95) # 95 upper 95th percent

max(coverage$sampsize_95)
# 53 is my cuttoff for the data

quantile(coverage$sampsize_95, 0.75)
##############################
# Apply threshold to the data

# load summary data
summary <- read.csv("5yr_summary/global_richness_summary.csv")

summary_filt <- summary %>% filter(number_checklists >= 53)
# went from 1.8 mil to 85,485

# now need to put this in raster form and merge with the urbanization raster
GHSLreproj<-rast("SMOD_global/reprojected.SMOD_global.tif")
#head(values(GHSLreproj))

# need to extract cell numbers, urbanization scores, and x and y coordinates from the raster
# turn into a data frame (this does not because the memory is not large enough)
# GHSL_df <- as.data.frame(GHSLreproj, xy=TRUE, cells=TRUE, na.rm=TRUE)


summary_filt$x <- xFromCell(GHSLreproj, summary_filt$cell) # extract the coordinates from the cells
summary_filt$y <- yFromCell(GHSLreproj, summary_filt$cell)

summary_filt$urban <- as.data.frame(terra::extract(GHSLreproj, summary_filt[,c(20:21)]))$SMOD_global

hist(summary_filt$x)

# definitely concentrated in the northern hemisphere
hist(summary_filt$y)
hist(summary_filt$urban)
summary_filt %>% group_by(urban) %>% summarise(n=n())

# remove ones with urbanization score of 10 (water)
summary_filt <- summary_filt %>% filter(!urban==10)
# save the thresholded data
#write.csv(summary_filt, "5yr_summary/summary_thresholded.csv", row.names=FALSE)




