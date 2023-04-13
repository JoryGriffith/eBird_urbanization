### Filtering for South American occurences

library(terra)
library(sf)
library(auk)
library("rnaturalearth")
library("rnaturalearthdata")
library(beepr)

# load global urbanization data
GHSL<-rast("SMOD_global/SMOD_global.tif")

# Load south american polygon
SAmerica <- ne_countries(scale = "medium", continent="south america", returnclass = "sf")
plot(SAmerica)
st_bbox(SAmerica)

eBird_2021 <- auk_ebd("eBird_2021_data/ebd_2021.txt")

###### filter ebird data to South America ######

eBird_2021 %>% 
  auk_bbox(SAmerica) %>% # filter out square
  auk_filter(select = cols, file = "eBird_2021_data/ebd_2021_SAm_unfilt.txt") # make dataframe

# then filter out columns I want
cols <- c("latitude", "longitude", "group identifier", "sampling event identifier",
          "scientific name", "observation count", "observer_id")

auk_ebd("eBird_2021_data/ebd_2021_SAm_unfilt.txt") %>% 
  auk_select(select = cols, file = "eBird_2021_data/ebd_2021_SAm_filt.txt")


##### Filtering for asian occurences
asia <- ne_countries(scale = "medium", continent="asia", returnclass = "sf")
plot(asia)
st_bbox(asia)

eBird_2021 %>% 
  auk_bbox(asia) %>% # filter out square
  auk_filter(select = cols, file = "eBird_2021_data/ebd_2021_asia_unfilt.txt") # make dataframe

# then filter out columns I want
cols <- c("latitude", "longitude", "group identifier", "sampling event identifier",
          "scientific name", "observation count", "observer_id")

auk_ebd("eBird_2021_data/ebd_2021_asia_unfilt.txt") %>% 
  auk_select(select = cols, file = "eBird_2021_data/ebd_2021_asia_filt.txt")

beep()


#### Filtering for Australia
australia <- ne_countries(scale = "medium", continent="oceania", returnclass = "sf")
plot(australia_reproj)
australia_reproj <- st_transform(australia, crs = 4326)
st_bbox(australia)

eBird_2021 %>% 
  auk_bbox(bbox=c(-179.9999999, -54.74922, 179.9999999, 18.80679)) %>% # filter out square
  auk_filter(select = cols, file = "eBird_2021_data/ebd_2021_australia_unfilt.txt") # make dataframe

# then filter out columns I want
cols <- c("latitude", "longitude", "group identifier", "sampling event identifier",
          "scientific name", "observation count", "observer_id")

auk_ebd("eBird_2021_data/ebd_2021_australia_unfilt.txt") %>% 
  auk_select(select = cols, file = "eBird_2021_data/ebd_2021_australia_filt.txt")

beep()








