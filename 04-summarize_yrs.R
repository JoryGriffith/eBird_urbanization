################### New data summary that accounts for repeat measures between years ####################
library(tidyverse)
library(lubridate)

######### Redoing the richness calculations over all of the years
years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

# loop for each square (skipped 5 because it's too big)
for (j in 8:length(names)){
datalist = vector("list", length = length(years))
  # loop for each year
  for (i in 1:length(years)) {
    dat <- read.table(paste("/Volumes/Backup/eBird/eBird_", years[i], "_data/custom_bbox/", names[j], "_", years[i], "_filt.txt", sep=""), 
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
months <- dat2 %>% group_by(cell, month) %>% 
  summarise(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)))

months_wide <- pivot_wider(months, id_cols = cell, names_from=month, values_from = number_checklists)
months_wide[is.na(months_wide)]<-0

## aggregate to get number of checklists and raw richness per cell
dat_summary <- dat2 %>%
  group_by(cell) %>%
  summarize(number_checklists=length(unique(SAMPLING.EVENT.IDENTIFIER)),
            total_SR=length(unique(SCIENTIFIC.NAME)),
            total_duration=sum(DURATION.MINUTES),
            avg_duration=mean(DURATION.MINUTES),
            no_months=length(unique(month))
  )

dat_summary2 <- merge(dat_summary, months_wide, by="cell")
dat_summary2$square=names[j]
# save as csv
write.csv(dat_summary2, paste("5yr_summary/", names[j], "_SR.csv", sep=""))
print(paste("finished", names[j]))

rm(dat)
rm(dat2)
rm(dat_summary)
rm(dat_summary2)
}




######### put all together and save
list_csv_files <- list.files(path = "5yr_summary/", pattern="*.csv")
#dat <- readr::read_csv(paste("5yr_summary/", list_csv_files, sep=""), id = "file_name")

names <- tolower(gsub('_SR.csv', "", list_csv_files))


for(i in 1:length(list_csv_files)) {                              # Head of for-loop
  assign(names[i],                                   # Read and store data frames
         read.csv(paste("5yr_summary/", list_csv_files[i], sep="")))
}

dat <- bind_rows(r1c1, r1c2, r1c3, r1c4, 
                 r2c1, r2c2aa, r2c2aba, r2c2abb, r2c2b, r2c3, r2c4, 
                 r3c1, r3c2, r3c3, r3c4, 
                 r4c1, r4c2, r4c3, r4c4)
#save all the summaries as a csv
write.csv(dat, "global_richness_summary.csv", row.names=FALSE)

# find top 500 cells
top_cells <- dat %>% slice_max(total_SR, n=500)
write.csv(top_cells, "top_500_cells.csv", row.names=FALSE)




