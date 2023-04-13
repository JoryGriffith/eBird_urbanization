########## In this script, I will take the cells with the highest richness and plot the species accumulation curves
### I want to get a good idea of how much sampling I need to get good curves in this high richness area

library(terra)
library(tidyverse)
library(sf)
library(auk)
library(iNEXT)
library(beepr)

GHSLreproj<-rast("SMOD_global/reprojected.SMOD_global.tif")

# Load eBird data
# only loading data with cells at the equator

##### R3C1 ######
# start with R3C1 because it doesnt have that much data
r3c1 <- read_ebd("eBird_2021_data/custom_bbox/R3C1_2021_filt.txt")
beep()
r3c1_vect <- vect(r3c1, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
GHSL_r3c1 <- crop(GHSLreproj, ext(r3c1_vect))
plot(GHSL_r3c1) # yes this is mostly water lol

# load summary data
r3c1_summary <- read.csv("eBird_2021_data/custom_bbox/summary/r3c1_summary.csv")
r3c1_summary2 <- r3c1_summary[which(r3c1_summary$number_checklists>1),]
# when you remove rows with only one checklist there is a definite loss of data

# this is not a very data rich area clearly, the highest richness is 48 and the highest checklists is 114

# the highest richness area is this little island off the coast of ecuador so that makes sense that it would be lower
# I think this block is mostly 

#### use the iNEXT package to look at coverage and plot the species accumulation curves

# extract cell numbers and values
xy=geom(r3c1_vect)
# get cell number that each point is in
r3c1$cell=cellFromXY(GHSL_r3c1, xy[,3:4])


temp <- r3c1 %>% 
  filter(cell=="39537153") %>% 
  group_by(scientific_name, sampling_event_identifier) %>% 
  summarize(present=n()) %>%
  mutate(present=1) %>%
  pivot_wider(names_from=sampling_event_identifier, values_from=present, values_fill=0) %>%
  column_to_rownames(var="scientific_name") %>%
  as.data.frame()

# convert this dataframe into data format for iNext
temp_inext <- as.incfreq(temp)

#and now using estimateD to get the species richness (q=0)
out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq") # none of this will work if there is only 1 checklist

ggiNEXT(out.inc, type=1) # plot of extrapolated species richness
ggiNEXT(out.inc, type=2) # plot of extrapolated sample coverage 
ggiNEXT(out.inc, type=3) # plot of species diversity and sample coverage

out.inc$iNextEst$coverage_based
out.inc$iNextEst$coverage_based %>% filter(Method=="Observed") 
# sample coverage for the cell with highest sp richness (1073362) is 0.96
# sample coverage for cell with most checklists (39537153) is 0.995
# cell with highest duration (15227359): 0.997
# cell with 17 checklists (39787660): 0.95
# cell with 6 checklists (39859213): 0.899
# cell with 4 checklists and SR of 17 (15817731): 0.92
# some of them are giving me a sample coverage of 1 for some reason
# should maybe try this with estimates of abundance as well instead of the incidence fgrequency thing??


###### R2C1 #######
r2c1 <- read_ebd("eBird_2021_data/custom_bbox/R2C1_2021_filt.txt")
beep()


r2c1_vect <- vect(r2c1, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
beep()
GHSL_r2c1 <- crop(GHSLreproj, ext(r2c1_vect))
plot(GHSL_r2c1) # southern part of the US and mexico

# load summary data
r2c1_summary <- read.csv("eBird_2021_data/custom_bbox/summary/r2c1_summary.csv")
r2c1_summary2 <- r2c1_summary[which(r2c1_summary$number_checklists>1),]
# when you remove rows with only one checklist there is a definite loss of data

# this is not a very data rich area clearly, the highest richness is 48 and the highest checklists is 114

# the highest richness area is this little island off the coast of ecuador so that makes sense that it would be lower
# I think this block is mostly 

#### use the iNEXT package to look at coverage and plot the species accumulation curves

# extract cell numbers and values
xy=geom(r2c1_vect)
# get cell number that each point is in
r2c1$cell=cellFromXY(GHSL_r2c1, xy[,3:4])

# random sample of 10 cells 
# (want to figure out how to randomly sample cells then calculate inc-freq data then put into separate lists to plot cverage side by side)
#cells <- sample(r2c1$cell, size=10, replace=FALSE)

temp <- r2c1 %>% 
  filter(cell == "812656") %>% # filter by some random grid cells for comparison
  group_by(scientific_name, sampling_event_identifier) %>% 
  summarize(present=n()) %>%
  mutate(present=1) %>%
  pivot_wider(names_from=sampling_event_identifier, values_from=present, values_fill=0) %>%
  column_to_rownames(var="scientific_name") %>%
  as.data.frame()

# convert this dataframe into data format for iNext
temp_inext <- as.incfreq(temp)

#and now using estimateD to get the species richness (q=0)
out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq") # none of this will work if there is only 1 checklist

ggiNEXT(out.inc, type=1, color.var="Order.q") # plot of extrapolated species richness
ggiNEXT(out.inc, type=2, color.var="Order.q") # plot of extrapolated sample coverage 
ggiNEXT(out.inc, type=3, color.var="Order.q") # plot of species diversity and sample coverage

#out.inc$iNextEst$coverage_based
out.inc$iNextEst$coverage_based %>% filter(Method=="Observed") 

# Sample coverage of cell w highest # checklists (812656): 0.999
# Sample coverage of highest species richness (16617236): 0.999
# sample coverage of random cell with only 2 checklists (5485): 0.67
# sample coverage of random cell with 4 checklists (5689): 0.54


########### Methods for creating a threshold of top cells #########
# load csv
#top_cells <- read.csv("500_richest_cells.csv")

# But now I need to load in the raw data and extract the actual bird sightings so 
# that I can look at species accum curves
unique(top_cells$square)
# r2c2AA, r2c2ABA, r2c3, r3c2, r2c1, r2c2ABB, r1c2

# but also need to do this for the different years uh oh this is going to be difficult

# maybe just start with 2021 as proof of concept?







##### Finding top squares for 2021 data
list_csv_files <- list.files(path = "eBird_2021_data/custom_bbox/summary", pattern="*.csv")
#
names <- tolower(gsub('_summary.csv', '_2021', list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("eBird_2021_data/custom_bbox/summary/", list_csv_files[i], sep="")))
}

df_2021<-bind_rows(r1c1_2021, r1c2_2021, r1c3_2021, r1c4_2021,
      r2c1_2021, r2c2aa_2021, r2c2aba_2021, r2c2abb_2021, r2c2b_2021, r2c3_2021, r2c4_2021,
      r3c1_2021, r3c2_2021, r3c3_2021, r3c4_2021,
      r4c1_2021, r4c2_2021, r4c3_2021, r4c4_2021)
df_2021<-df_2021[c(2:23)] # take away random X variable thing

# remove to save space
remove(r1c1_2021, r1c2_2021, r1c3_2021, r1c4_2021,
       r2c1_2021, r2c2aa_2021, r2c2aba_2021, r2c2abb_2021, r2c2b_2021, r2c3_2021, r2c4_2021,
       r3c1_2021, r3c2_2021, r3c3_2021, r3c4_2021,
       r4c1_2021, r4c2_2021, r4c3_2021, r4c4_2021)

top_cells_2021 <- df_2021 %>% slice_max(total_SR, n=500)
# find square with most top cells and start with that
top_cells_2021 %>% group_by(square) %>% summarise(n=n())
# r2c2AA is the one with the most top cells
# now need to load it and filter out the squares that are at the top

#### Load r2c2AA
# filter top cells by r2c2AA
top_cells_r2c2AA <- top_cells_2021 %>% filter(square=="r2c2AA")

dat <- read_ebd("eBird_2021_data/custom_bbox/R2c2AA_2021_filt.txt")
beep()
vect <- vect(dat, crs=crs(GHSLreproj),geom=c("longitude","latitude"))
GHSL_crop <- crop(GHSLreproj, ext(vect))

xy=geom(vect)

# get cell number that each point is in
dat$cell=cellFromXY(GHSL_crop, xy[,3:4])

# filter by cell numbers that I want
dat_filt <- dat %>% filter(cell %in% top_cells_r2c2AA$cell)


##### calculate sample completeness for all
# split by cell
ebird.split <- dat_filt %>% group_by(cell) %>% group_split()

# make data frame for output
output <- data.frame(matrix(NA,
                            nrow=length(ebird.split),
                            ncol=9))
# make column names for output
names(output)<-c("cell_ID", "SC", "obs.richness", "obsLCL", "obsUCL", "method", "st.richness", "st.LCL", "st.UCL")

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
    out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq")
    out.inc.filt<-out.inc$iNextEst$coverage_based %>% filter(Method=="Observed")  # filter for only observed method (for sample size)
    
    sampling_profile <- estimateD(temp_inext, q=0, datatype="incidence_freq",base="coverage", 
                                  level=0.80, conf=0.95) # set to a coverage level of 60
    
    # add output to a data frame
    output$SC[i] <- out.inc.filt$SC # sample coverage
    output$obs.richness[i] <- out.inc.filt$qD
    output$obsLCL[i] <- out.inc.filt$qD.LCL
    output$obsUCL[i] <- out.inc.filt$qD.UCL
    output$method[i] <- sampling_profile$Method
    output$st.richness[i] <- sampling_profile$qD
    output$st.LCL[i] <- sampling_profile$qD.LCL
    output$st.UCL[i] <- sampling_profile$qD.UCL
    # add new columns to data frame
  }
},
warning = function(w){
  output$warnings[i] <<- w$message
  invokeRestart("muffleWarning")
})


### Try sampling one checklist at a time until I read 95%

# to do so, I will sample one checklist at a time and calculate coverage and stop when I reach 95%

# first need to turn into incidence frequency data

write.csv(output, "thresholding/r2c2AA_2021_coverage.csv")


# now, I need to write some crazy code to do this bootstrapping approach


# Choose a random square and run iNEXT on it
temp <- dat_filt %>% 
  filter(cell==1517104) %>% # select list element
  group_by(scientific_name, sampling_event_identifier) %>% 
  summarize(present=n()) %>%
  mutate(present=1) %>%
  pivot_wider(names_from=sampling_event_identifier, values_from=present, values_fill=0) %>%
  column_to_rownames(var="scientific_name") %>%
  as.data.frame() 

# convert this dataframe into data format for iNext
temp_inext <- as.incfreq(temp)
#and now using estimateD to get qD
out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=200, nboot=100)
?iNEXT
out.inc$iNextEst



####################
# Some summary stats of the top cells
top_cells <- read.csv("500_richest_cells.csv")
GHSLreproj<-rast("SMOD_global/reprojected.SMOD_global.tif")
plot(GHSLreproj)

ggplot(top_cells, aes(x=x, y=y))+
  geom_point()
# they are all nearish the equator, that's good

# try to rasterize
top_cells_vect <- vect(top_cells, crs=crs(GHSLreproj), geom=c("x","y"))
summary(values(top_cells_vect))
richness.raster <- terra::rasterize(top_cells_vect, GHSLreproj, field="total_SR", fun = mean)

plot(richness.raster)
values(richness.raster)
# don't know why this is not working, come back to it


###START HERE






#########################################
# Need to extract all of the raw data from the squares from all of the years and put together
# So that I can then calculate species accumulation curves for it
# I will do it by square

# for some reason it is 501 so let's cut off the last one and the column at the beginning
top_cells <- top_cells[c(1:500), c(2:23)]
unique(top_cells$square)
# "r2c1" "r2c2AA" "r2c2ABA" "r2c2ABB"  "r2c2B" "r2c3"   "r2c4"  "r3c2"  "r3c3"  "r3c4"  
# Load each one of these separately in a loop

top_cells %>% group_by(square) %>% summarise(n=n())
# r3c4 has the least number of squares so will start with that

years <- c(2017, 2018, 2019, 2020, 2021)


##### r2c1
top_cells_r2c1 <- top_cells %>% filter(square=="r2c1")
# make list
datalist_r2c1 = vector("list", length = length(years))

for (i in 1:5) {
  dat <- read.table(paste("eBird_", years[i], "_data/custom_bbox/r2c1_", years[i], "_filt.txt", sep=""), header=TRUE)
  
  dat_filt <- dat %>% filter(cell %in% top_cells_r2c1$cell)
  
  datalist_r2c1[[i]] <- dat_filt
  #write.csv(dat_filt, paste("thresholding/r2c1_", years[i], "_topcells.csv", sep=""))
}
# bind lists together
r2c1_data <- dplyr::bind_rows(datalist_r2c1)
# save as csv
write.csv(r2c1_data, "thresholding/r2c1_topcells.csv")




##### R3C4
top_cells_r3c4 <- top_cells %>% filter(square=="r3c4")
# make list
datalist_r3c4 = vector("list", length = length(years))

for (i in 1:5) {
  dat <- read.table(paste("eBird_", years[i], "_data/custom_bbox/r3c4_", years[i], "_filt.txt", sep=""), header=TRUE)
  
  dat_filt <- dat %>% filter(cell %in% top_cells_r3c4$cell)
  
  datalist_r3c4[[i]] <- dat_filt
  #write.csv(dat_filt, paste("thresholding/r3c4_", years[i], "_topcells.csv", sep=""))
}
# bind lists together
r3c4_data <- dplyr::bind_rows(datalist_r3c4)
# save as csv
write.csv(r3c4_data, "thresholding/r3c4_topcells.csv")






##########################
# Write loop to calculate the coverage and find the number of checklists that is closest to 95% coverage

# load data
#squares <- c("r2c2AA", "r2c2ABA", "r2c3", "r3c2", "r2c1", "r2c2ABB" ,"r1c2")

# make datalist for output
datalist = vector("list", length = length(squares))

for (j in 7:length(squares)) {

dat_top <- read.csv(paste("thresholding/", squares[j], "_topcells.csv", sep=""))

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
    out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=500, nboot=100, endpoint = 500)
    out.inc.obs <- out.inc$DataInfo
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
  }
},
warning = function(w){
  output$warnings[i] <<- w$message
  invokeRestart("muffleWarning")
})

output$square <- square[j]
datalist[[j]] <- output

}








############
# Test

dat_top <- read.csv("thresholding/r3c4_topcells.csv")

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
    out.inc <- iNEXT(temp_inext, q=0, datatype="incidence_freq", knots=100, nboot=50, endpoint = 500)
    out.inc.obs <- out.inc$DataInfo
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
  }
},
warning = function(w){
  output$warnings[i] <<- w$message
  invokeRestart("muffleWarning")
})







