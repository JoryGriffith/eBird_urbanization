## This is a script to filter the data by seasons so I can conduct seasonal analyses
library(auk)

cols <- c("latitude", "longitude", "group identifier", "sampling event identifier",
          "scientific name", "observation count", "observer_id", "observation_date", "duration_minutes")

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

# Summer
for (j in 1:length(years)){
  for (i in 1:length(names)){
    # filter for dates I want
    auk_ebd(paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/custom_bbox/", names[i], "_", years[j], "_unfilt.txt", sep="")) %>% 
      auk_date() %>% # add summer date to filter for 
      auk_filter(file=paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_unfilt.txt", sep=""))
    
    # then filter for columns I want
    auk_ebd(paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_unfilt.txt", sep="")) %>% 
      auk_select(file=paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/summer/", names[i], "_", years[j], "_summer_filt.txt", sep=""))
  }
}


# Winter
for (j in 1:length(years)){
  for (i in 1:length(names)){
    # filter for dates I want
    auk_ebd(paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/custom_bbox/", names[i], "_", years[j], "_unfilt.txt", sep="")) %>% 
      auk_date() %>% # add winter date to filter for 
      auk_filter(file=paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_unfilt.txt", sep=""))
    
    # then filter for columns I want
    auk_ebd(paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_unfilt.txt", sep="")) %>% 
      auk_select(file=paste("/Volumes/Expansion/eBird/eBird_", years[j], "_data/winter/", names[i], "_", years[j], "_winter_filt.txt", sep=""))
  }
}






