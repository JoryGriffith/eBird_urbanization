##### Combining GHM and GHSL layers
library(terra)
library(tidyverse)

# Load GHSL layer
GHSLreproj<-rast("SMOD_global/reprojected.SMOD_global.tif")

# Load GHM layer
#GHM <- rast("gHM/gHM.tif")
## Reproject 
#GHMreproject <- project(GHM, y="+proj=longlat +datum=WGS84", method="near")
## save as raster
#writeRaster(GHMreproject, "ghM/GHMreprojected.tif", overwrite=TRUE)

# load GHM
GHM <- rast("ghM/GHMreprojected.tif")

# assign values more than 0.5 an NA 
#GHM[GHM$gHM>=0.5] <- NA

crs(GHM)==crs(GHSLreproj)
# they are in the same crs
# just need to make them line up

# try with nearest neighbor interpolation
dat <- resample(GHM, GHSLreproj, method="cubicspline")
# it should be aligned now and I can create a raster stack
# should look more into interpolation methods for estimating the new cell values in the GHM
plot(dat)

# make raster stack
GHSL

extent <- c(-100, -70, 30, 50)
colors=c("palegreen", "darkorange", "firebrick2")
GHSL.test <- GHSLreproj
GHSL.test[(dat > 0.5 & GHSL.test == 11)] <- NA # turn everything that is wilderness and over 0.5 as NA
GHSL.test[(GHSL.test ==10),] <- NA # turn water into NA
# need to refigure out this code because it doesn't seem to be working
plot(GHSL.test, ext=extent, col=colors)


# making sure it worked
#dat[(dat>0.5 & GHSL == 11)] <- 1000

#plot(dat)
#plot(GHSL)
#plot(GHSL.test)
# looks good

# change values to 1, 2, and 3
#GHSL.test2 <- GHSL.test
#GHSL.test2 <- subst(GHSL.test2, 10, NA) # turn water into NA
#GHSL.test2 <- subst(GHSL.test2, 10, NA) # turn wilderness into 1
#GHSL.test2 <- subst(GHSL.test2, 10, NA) # turn peri-urban into 2
#GHSL.test2 <- subst(GHSL.test2, 10, NA) # turn urban into 3

# make raster stack for saving
#urban_raster <- c(GHSL.test, dat)
# reproject

# save as tif file
writeRaster(GHSL.test, "SMOD_global/GHSL_filtered.tif")


