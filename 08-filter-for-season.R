## This is a script to filter the data by seasons so I can conduct seasonal analyses
library(auk)

cols <- c("latitude", "longitude", "group identifier", "sampling event identifier",
          "scientific name", "observation count", "observer_id", "observation_date", "duration_minutes", "effort_distance_km")

namesN <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4")

namesS <- c("r3c1", "r3c2", "r3c3", "r3c4",
            "r4c1", "r4c2", "r4c3", "r4c4")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

library(doParallel)
numCores<-detectCores()
cl <- makeCluster(numCores)
registerDoParallel(numCores)

# Summer N hemisphere
for (j in 1:length(years)){
  foreach (i=1:length(namesN)) %dopar%{
 # for (i in 1:length(namesN)){
    # filter for dates I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", namesN[i], "_", years[j], "_unfilt.txt", sep="")) %>% 
      auk_date(c("*-06-01", "*-08-31")) %>% # add summer date to filter for 
      auk_filter(file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", namesN[i], "_", years[j], "_summer_unfilt.txt", sep=""), overwrite=TRUE)
    
    # then filter for columns I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", namesN[i], "_", years[j], "_summer_unfilt.txt", sep="")) %>% 
      auk_select(select=cols, file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", namesN[i], "_", years[j], "_summer_filt.txt", sep=""), overwrite=TRUE)
    print(paste("finished", namesN[i]))
  }
print(paste("finished", years[j]))
  }


# Winter N hemisphere
for (j in 1:length(years)){
  foreach (i=1:length(namesN)) %dopar%{
    # filter for dates I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", namesN[i], "_", years[j], "_unfilt.txt", sep="")) %>% 
      auk_date(c("*-12-01", "*-02-29")) %>% # add winter date to filter for 
      auk_filter(file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", namesN[i], "_", years[j], "_winter_unfilt.txt", sep=""), overwrite=TRUE)
    
    # then filter for columns I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", namesN[i], "_", years[j], "_winter_unfilt.txt", sep="")) %>% 
      auk_select(select=cols, file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", namesN[i], "_", years[j], "_winter_filt.txt", sep=""), overwrite=TRUE)
    print(paste("finished", namesN[i]))
  }
  print(paste("finished", years[j]))
}


## Winter S hemisphere
for (j in 1:length(years)){
  foreach (i=1:length(namesS)) %dopar%{
#  for (i in 1:length(namesS)){
    # filter for dates I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", namesS[i], "_", years[j], "_unfilt.txt", sep="")) %>% 
      auk_date(c("*-06-01", "*-08-31")) %>% # add summer date to filter for 
      auk_filter(file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", namesS[i], "_", years[j], "_winter_unfilt.txt", sep=""), overwrite=TRUE)
    
    # then filter for columns I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", namesS[i], "_", years[j], "_winter_unfilt.txt", sep="")) %>% 
      auk_select(select=cols, file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/winter/", namesS[i], "_", years[j], "_winter_filt.txt", sep=""), overwrite=TRUE)
    print(paste("finished", namesS[i]))
  }
  print(paste("finished", years[j]))
}

## Summer S hemisphere
for (j in 1:length(years)){
  foreach (i=1:length(namesS)) %dopar%{
    # filter for dates I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/custom_bbox/", namesS[i], "_", years[j], "_unfilt.txt", sep="")) %>% 
      auk_date(c("*-12-01", "*-02-29")) %>% # add winter date to filter for 
      auk_filter(file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", namesS[i], "_", years[j], "_summer_unfilt.txt", sep=""), overwrite=TRUE)
    
    # then filter for columns I want
    auk_ebd(paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", namesS[i], "_", years[j], "_summer_unfilt.txt", sep="")) %>% 
      auk_select(select=cols, file=paste("/Volumes/Backup/eBird/eBird_", years[j], "_data/summer/", namesS[i], "_", years[j], "_summer_filt.txt", sep=""), overwrite=TRUE)
    print(paste("finished", namesS[i]))
  }
  print(paste("finished", years[j]))
}

