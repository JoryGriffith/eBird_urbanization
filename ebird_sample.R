##### Practicing with the eBird sample data ########
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)


# load sample data
ebd_sample<-read_ebd("ebd_US-AL-101_202204_202204_relApr-2022_SAMPLE/ebd_US-AL-101_202204_202204_relApr-2022.txt")
nrow(ebd_sample)
# 1339

# exploratory analysis and visualization

ebird_sf <- ebd_sample %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(scientific_name)

# plot
plot(ebird_sf)
# all of these are from alabama

# download the map of alabama
library(maptools)
library(mapdata)

states <- map_data("state")
AL <- states %>% subset(region=="alabama")

ggplot() + 
 # geom_polygon(data = AL, aes(x=long, y = lat, group = group),
  #             color= "black", fill= "white")+
  geom_sf(data = ebird_sf, size=0.2)+
  coord_sf()
# they are all in a very small part of alabama
# this is not going to be super helpful to work with

ggplot() + 
  geom_polygon(data = Japan, aes(x=long, y = lat, group = group),
               color= "black", fill= "white") +
  geom_sf(data = jpn, size=0.2, color = "red", fill = "white") + 
  coord_sf()








