####### This is a script to look at the coverage of the aggregated 2021 data fo the globe #####
########## NOT USING ANYMORE BECAUSE I MADE A MISTAKE IN THESE SUMMARIES ######
library(automap)
library(nlme)
library(spdep)
library(ncf)
library(tidyverse)
library(readr)
#### Took all the data from 2021 and the aggregated it by grid cell 
  # and calculated the species richness and number of checklists for each grid cell
# had to do it in chunks because there was so much data
# but now can load it all and combine it to take a better look


###### Load all 2021 Data
#list_csv_files <- list.files(path = "eBird_2021_data/custom_bbox/summary", pattern="*.csv")
#
#names <- tolower(gsub('_summary.csv', '_2021', list_csv_files))
#
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         read.csv(paste("eBird_2021_data/custom_bbox/summary/", list_csv_files[i], sep="")))
#}
#
#df_2021<-bind_rows(r1c1_2021, r1c2_2021, r1c3_2021, r1c4_2021,
#      r2c1_2021, r2c2aa_2021, r2c2aba_2021, r2c2abb_2021, r2c2b_2021, r2c3_2021, r2c4_2021,
#      r3c1_2021, r3c2_2021, r3c3_2021, r3c4_2021,
#      r4c1_2021, r4c2_2021, r4c3_2021, r4c4_2021)
#df_2021<-df_2021[c(2:23)] # take away random X variable thing
#
## remove to save space
#remove(r1c1_2021, r1c2_2021, r1c3_2021, r1c4_2021,
#       r2c1_2021, r2c2aa_2021, r2c2aba_2021, r2c2abb_2021, r2c2b_2021, r2c3_2021, r2c4_2021,
#       r3c1_2021, r3c2_2021, r3c3_2021, r3c4_2021,
#       r4c1_2021, r4c2_2021, r4c3_2021, r4c4_2021)
#
#
###### Load 2020 data
#list_csv_files <- list.files(path = "eBird_2020_data/custom_bbox/summary", pattern="*.csv")
#
#names <- gsub('_summary.csv', '_2020', list_csv_files)
#
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         read.csv(paste("eBird_2020_data/custom_bbox/summary/", list_csv_files[i], sep="")))
#}
#
#df_2020<-bind_rows(r1c1_2020, r1c2_2020, r1c3_2020, r1c4_2020,
#                   r2c1_2020, r2c2AA_2020, r2c2ABA_2020, r2c2ABB_2020, r2c2B_2020, r2c3_2020, r2c4_2020,
#                   r3c1_2020, r3c2_2020, r3c3_2020, r3c4_2020,
#                   r4c2_2020, r4c3_2020)
#df_2020<-df_2020[c(2:23)] # take away random X variable thing
#
## remove to save space
#remove(r1c1_2020, r1c2_2020, r1c3_2020, r1c4_2020,
#       r2c1_2020, r2c2AA_2020, r2c2ABA_2020, r2c2ABB_2020, r2c2B_2020, r2c3_2020, r2c4_2020,
#       r3c1_2020, r3c2_2020, r3c3_2020, r3c4_2020,
#       r4c2_2020, r4c3_2020, r4c4_2020)
#
#
#
#
##### Load 2019 data
#list_csv_files <- list.files(path = "eBird_2019_data/custom_bbox/summary", pattern="*.csv")
#
#names <- gsub('_summary.csv', '_2019', list_csv_files)
#
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         read.csv(paste("eBird_2019_data/custom_bbox/summary/", list_csv_files[i], sep="")))
#}
#
#df_2019<-bind_rows(r1c1_2019, r1c2_2019, r1c3_2019, r1c4_2019,
#                   r2c1_2019, r2c2AA_2019, r2c2ABA_2019, r2c2ABB_2019, r2c2B_2019, r2c3_2019, r2c4_2019,
#                   r3c1_2019, r3c2_2019, r3c3_2019, r3c4_2019,
#                   r4c2_2019, r4c4_2019)
#df_2019<-df_2019[c(2:23)] # take away random X variable thing
#
## remove to save space
#remove(r1c1_2019, r1c2_2019, r1c3_2019, r1c4_2019,
#       r2c1_2019, r2c2AA_2019, r2c2ABA_2019, r2c2ABB_2019, r2c2B_2019, r2c3_2019, r2c4_2019,
#       r3c1_2019, r3c2_2019, r3c3_2019, r3c4_2019,
#       r4c2_2019, r4c4_2019)
#
#
##### Load 2018
#list_csv_files <- list.files(path = "eBird_2018_data/custom_bbox/summary", pattern="*.csv")
#
#names <- gsub('_summary.csv', '_2018', list_csv_files)
#
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         read.csv(paste("eBird_2018_data/custom_bbox/summary/", list_csv_files[i], sep="")))
#}
#
#df_2018<-bind_rows(r1c1_2018, r1c2_2018, r1c3_2018, r1c4_2018,
#                   r2c1_2018, r2c2AA_2018, r2c2ABA_2018, r2c2ABB_2018, r2c2B_2018, r2c3_2018, r2c4_2018,
#                   r3c1_2018, r3c2_2018, r3c3_2018, r3c4_2018,
#                   r4c2_2018, r4c3_2018, r4c4_2018)
#df_2018<-df_2018[c(2:23)] # take away random X variable thing
#
## remove to save space
#remove(r1c1_2018, r1c2_2018, r1c3_2018, r1c4_2018,
#       r2c1_2018, r2c2AA_2018, r2c2ABA_2018, r2c2ABB_2018, r2c2B_2018, r2c3_2018, r2c4_2018,
#       r3c1_2018, r3c2_2018, r3c3_2018, r3c4_2018,
#       r4c2_2018, r4c3_2018, r4c4_2018)
#
#
###### Load 2017
#list_csv_files <- list.files(path = "eBird_2017_data/custom_bbox/summary", pattern="*.csv")
#
#names <- gsub('_summary.csv', '_2017', list_csv_files)
#
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         read.csv(paste("eBird_2017_data/custom_bbox/summary/", list_csv_files[i], sep="")))
#}
#
#df_2017<-bind_rows(r1c1_2017, r1c2_2017, r1c3_2017, r1c4_2017,
#                   r2c1_2017, r2c2AA_2017, r2c2ABA_2017, r2c2ABB_2017, r2c2B_2017, r2c3_2017, r2c4_2017,
#                   r3c1_2017, r3c2_2017, r3c3_2017, r3c4_2017,
#                   r4c2_2017,  r4c4_2017)
#df_2017<-df_2017[c(2:23)] # take away random X variable thing
#
## remove to save space
#remove(r1c1_2017, r1c2_2017, r1c3_2017, r1c4_2017,
#       r2c1_2017, r2c2AA_2017, r2c2ABA_2017, r2c2ABB_2017, r2c2B_2017, r2c3_2017, r2c4_2017,
#       r3c1_2017, r3c2_2017, r3c3_2017, r3c4_2017,
#       r4c2_2017, r4c4_2017)
#
#
### Now I have all the data frames loaded and I want to turn them into a large data frame by summarizing by raster square
#df_2017$year <- 2017
#df_2018$year <- 2018
#df_2019$year <- 2019
#df_2020$year <- 2020
#df_2021$year <- 2021
#
## bind together
#df <- bind_rows(df_2017, df_2018, df_2019, df_2020, df_2021) #3506692 
#
## summarise by raster cell and square
#
#df_sum <- df %>% group_by(cell, square) %>% summarise(x=mean(x),
#                                                      y=mean(y),
#                                                      cell=mean(cell),
#                                                      SMOD_global=mean(SMOD_global),
#                                                      number_checklists = sum(number_checklists),
#                                                      total_SR = sum(total_SR),
#                                                      total_duration = sum(total_duration),
#                                                      avg_duration = mean(avg_duration),
#                                                      no_months=sum(no_months),
#                                                      square=unique(square),
#                                                      jan = sum(X1), # number of checklists in each month
#                                                      feb=sum(X2),
#                                                      mar=sum(X3),
#                                                      apr=sum(X4),
#                                                      may=sum(X5),
#                                                      jun=sum(X6),
#                                                      jul=sum(X7),
#                                                      aug=sum(X8),
#                                                      sep = sum(X9),
#                                                      oct = sum(X10),
#                                                      nov=sum(X11),
#                                                      dec=sum(X12),
#                                                      )
#
#
#
#
#
#
## save as csv
#write.csv(df_sum, "summary_allyears.csv", row.names=F)
#

######### Load data frame
df <- read.csv("summary_allyears.csv")

# filter for 5 checklists
df_5 <- df[which(df$number_checklists>5), ]



summary(df_5)
# most are in the northern hemisphere, not as many near the equator
hist(df_5$y)

hist(df_5$SMOD_global)

hist(df_5$total_SR)

# filter for 10 checklists
df_10 <- df[which(df$number_checklists>10), ]

# there are 2,186,146 squares with checklists from 2017-2021
summary(df)
# There is a good spread of longitude, wit points from 179 to -179
# there is also a good spread of latitude, with points from 80 to -81
# species richness goes from 1 to 373


# spread of latitude 
hist(df$y)
# most of the data is in the northern hemisphere and most is at more northern latitudes

# spread of longitude
hist(df$x)

# spread of urbanization
hist(df$SMOD_global)
# mostly non-urbanized

# spread of number of checklists
# mostly small number
hist(log(df$number_checklists))


hist(df$total_SR)
# very log normal
hist(log(df$total_SR))



# relationship between number of checklists and richness
ggplot(df, aes(x=log(number_checklists), y=log(total_SR)))+
  geom_point() + 
  geom_smooth()+
  labs(x= "log # checklists", 
       y= "log total SR")
# definitely a positive relationship
# asymptotic as richness asymptotes



# relationship between richness and latitude
ggplot(df, aes(x=y, y=total_SR))+
  geom_point() + 
  geom_smooth()+
  labs(x= "latitude", 
       y= "total SR")
# nothing super strong here

# greater than 5 checklists
ggplot(df_5, aes(x=y, y=total_SR))+
  geom_point() + 
  geom_smooth()+
  labs(x= "latitude", 
       y= "total SR")

# more than 10 checklists
ggplot(df, aes(x=y, y=log(total_SR)))+
  geom_point() + 
  geom_smooth()+
  labs(x= "latitude", 
       y= "total SR")

# urbanization score and log richness
ggplot(df, aes(x=SMOD_global, y=log(total_SR)))+
  geom_point() + 
  geom_smooth(method="lm")+
  labs(x= "urbanization", 
       y= "log SR")
# positive relationship


# urbanization score and checklists
ggplot(df, aes(x=SMOD_global, y=log(number_checklists)))+
  geom_point() + 
  geom_smooth(method="lm")+
  labs(x= "urbanization", 
       y= "log # checklists")
# positive relationship


# latitude and number of checklists
ggplot(df, aes(x=y, y=log(number_checklists)))+
  geom_point() + 
  geom_smooth()+
  labs(x= "latitude", 
       y= "log # checklists")
# not really a pattern there, that's good


# latitude and richness by urbanization score
ggplot(df, aes(x=y, y=log(number_checklists)))+
  geom_point() + 
  geom_smooth()+
  facet_wrap(~SMOD_global)+
  labs(x= "latitude", 
       y= "log # checklists")


### Duration
# relationship between total duration and richness
ggplot(df, aes(x=log(total_duration), y=log(total_SR)))+
  geom_point() + 
  geom_smooth()+
  labs(x= "log duration", 
       y= "log total SR") # very strong relationship here

# urbanization score and duration
ggplot(df, aes(y=log(total_duration), x=log(SMOD_global)))+
  geom_point() + 
  geom_smooth(method="lm")+
  labs(y= "log duration", 
       x= "urbanization") # positive relationship

# total duration and number of checklists
ggplot(df, aes(y=log(total_duration), x=log(number_checklists)))+
  geom_point() + 
  geom_smooth(method="lm")+
  labs(y= "log duration", 
       x= "log number checklists") # very positive (obviously)


##### Month stats
# change all NAs to 0s
df[is.na(df)]<-0
# look at number of months
hist(df$no_months)
# most have only one month and the rest are pretty even

month_summary <- df %>% pivot_longer(cols=c(11:22))

month_summary %>% group_by(name) %>% summarise(sum=sum(value)) %>% 
  ggplot(aes(x=name, y=sum))+
  geom_col()
# relatively even spread but definitely the most checklists in april

# by square
month_summary %>% group_by(name, square) %>% filter(square=="r4c2") %>% summarise(sum=sum(value)) %>% 
  ggplot(aes(x=name, y=sum))+
  geom_col()+
  facet_wrap(~square)
# similar pattern across different areas





### Look at squares that are near the equator to see the sampling effort and richness
df_equator <- df %>% filter(y<1 & y>-1)
hist(df_equator$total_SR)
mean(df_equator$total_SR) # mean is 33
median(df_equator$total_SR) # median is 19
plot(log(df_equator$number_checklists)~log(df_equator$total_SR)) # very steep slope, makes sense

###### Run some linear models

# simple linear model with urbanization as a response
mod1 <- gls(log(total_SR) ~ y*SMOD_global + number_checklists + x, df_10)
plot(mod1) # need to also be thinking about the kind of model I want to do because its not a linear relationship
summary(mod1)

# check for spatial autocorrelation
df_10$mod1resids <- residuals(mod1, type="normalized")

residsI <- spline.correlog(x=df_10$x, y=df_10$y,
                           z=df_10$mod1Resids, resamp=50, quiet=TRUE)

# this doesn't seem to be working
datSPDF <- df_10
coordinates(datSPDF) <- c("x","y")
plot(autofitVariogram(mod1resids~1, input_data = datSPDF))
# ok this is working now but it still looks a bit weird, come back to this



####### Filter out 500 cells with highest richness for the thresholding

top_cells <- df %>% slice_max(total_SR, n=500)
# for some reason it is giving me 504 instead of 500

ggplot(top_cells, aes(x=x, y=y))+
  geom_point()
# a lot are at higher latitudes, probably because the sampling effort is so high there

write.csv(top_cells, file="500_richest_cells.csv")









