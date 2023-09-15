############## Calculating threshold season ############
library(terra)
library(tidyverse)
library(sf)
library(iNEXT)

# here I will calculate the threshold I need seperately for summer and winter months using the iNEXT package

######## Summer ######

summer <- read.csv("summer_richness_summary.csv")
summer_top_cells <- summer %>% slice_max(total_SR, n=500) # find top 500 cells
unique(summer_top_cells$square)

names <- c("r1c1", "r1c2", "r1c3", "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4", "r3c2", "r3c3", "r3c4") # figure out names

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# filter out cells that I want
for (j in 1:length(names)) {
  
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

names <- c("r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4", "r3c2", "r3c3") # figure out names

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# filter out cells that I want
for (j in 1:length(names)) {
  
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

for (j in 1:length(names)) {
  
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

