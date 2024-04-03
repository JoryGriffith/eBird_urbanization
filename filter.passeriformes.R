# Script for filtering for passerines
# starting from the summarise_years step
library(awk)
library(tidyverse)
library(lubridate)

years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")
# get taxonomy
tax <- get_ebird_taxonomy()
tax2 <- tax[,c(1,6)]


######### Redoing the richness calculations over all of the years
years <- c(2017, 2018, 2019, 2020, 2021, 2022)

names <- c("r1c1", "r1c2", "r1c3", "r1c4",
           "r2c1", "r2c2AA", "r2c2ABA", "r2c2ABB", "r2c2B", "r2c3", "r2c4",
           "r3c1", "r3c2", "r3c3", "r3c4",
           "r4c1", "r4c2", "r4c3", "r4c4")

j=18
i=4
# no passiformes in r4c1 or r4c3
# loop for each square (skipped 5 because it's too big)
for (j in 19:length(names)){
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
    dat.merged <- merge(dat, tax2, by.x="SCIENTIFIC.NAME", by.y="scientific_name") %>% 
      filter(order=="Passeriformes") # filter for passeriformes
    datalist[[i]] <- dat.merged # put in a list
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
write.csv(dat_summary2, paste("5yr_summary/", names[j], "_SR.passeriformes.csv", sep=""))
print(paste("finished", names[j]))

rm(dat)
rm(dat2)
rm(dat_summary)
rm(dat_summary2)
}
names
