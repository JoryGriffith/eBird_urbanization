#### Filtering different years of the eBird data
library(auk)
library(tidyverse)
library(sf)

# make list of years to include

# load whole dataset
ebd <- auk_ebd("ebd_relOct-2022/ebd_relOct-2022.txt",
               file_sampling = "ebd_sampling_relOct-2022/ebd_sampling_relOct-2022.txt")

years <- c(2017, 2018, 2019, 2020)
output_files <- c("eBird_2017_data/ebd_2017.txt", "eBird_2018_data/ebd_2018.txt", "eBird_2019_data/ebd_2019.txt", "eBird_2020_data/ebd_2020.txt")
output_sampling <- c("eBird_2017_data/ebd_sampling_2017.txt", "eBird_2018_data/ebd_sampling_2018.txt", "eBird_2019_data/ebd_sampling_2019.txt", "eBird_2020_data/ebd_sampling_2020.txt")


for(i in 1:1) {

  ebd %>% 
  auk_year(years[i]) %>% # filter for last year
  auk_complete() %>% # filter for complete checklists
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% # filter out incidental
  auk_distance(distance = c(0, 5)) %>%  # filter for distances less than 5km 
  auk_duration(duration = c(0, 300)) %>% # filter to less than 5 hours
  auk_filter(file = output_files[i], file_sampling = output_sampling[i], overwrite=TRUE) # create an output file from this

}

