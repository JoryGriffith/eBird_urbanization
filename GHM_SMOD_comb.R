##### Combining GHM and GHSL layers
library(terra)
library(tidyverse)

# Load GHSL layer
GHSL<-rast("/Volumes/Backup/eBird/SMOD_global/SMOD_global.tif")
plot(GHSL)
# Load GHM layer
GHM <- rast("/Volumes/Backup/eBird/gHM/gHM.tif")
plot(GHM)

# assign values more than 0.5 an NA 
#GHM[GHM$gHM>=0.5] <- NA

crs(GHM)==crs(GHSL)
# they are already in the same crs!
# just need to make them line up

# Make cells align with cubic spline interpolation
dat <- resample(GHM, GHSL, method="cubicspline")

# it should be aligned now and I can create a raster stack
# should look more into interpolation methods for estimating the new cell values in the GHM
plot(dat)


extent <- c(-100, -70, 30, 50)
colors=c("palegreen", "yellow", "pink", "darkorange", "blue", "firebrick2", "purple", "darkgreen")
GHSL.test <- GHSL
plot(GHSL.test)
GHSL.test[(dat > 0.5 & GHSL.test %in% c(11,12,13))] <- NA # turn everything that is natural and over 0.5 as NA
summary(values(GHSL.test))
GHSL.test[(GHSL.test==10)] <- NA # turn water into NA
GHSL[(GHSL==10)] <- NA

plot(GHSL, col=colors)
plot(GHSL.test, col=colors) # compare the two to make sure it worked, yes looks good!


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
writeRaster(GHSL.test, "/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif", overwrite=TRUE)

GHSL.test2 <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filtMollweide.tif")
plot(GHSL.test2)
GHSL.test2 <- subst(GHSL.test2, 10, NA) # turn water into NA
GHSL.test2 <- subst(GHSL.test2, 30, 3) # turn urban into 3
GHSL.test2 <- subst(GHSL.test2, 11, 1)
GHSL.test2 <- subst(GHSL.test2, 12, 1) 
GHSL.test2 <- subst(GHSL.test2, 13, 1) # turn 11,12, and 13 into natural
GHSL.test2 <- subst(GHSL.test2, 21, 2)
GHSL.test2 <- subst(GHSL.test2, 22, 2)
GHSL.test2 <- subst(GHSL.test2, 23, 2)

plot(GHSL.test2)

writeRaster(GHSL.test2, "/Volumes/Expansion/eBird/SMOD_global/GHSL_filt_3cat.tif")

filt <- rast("/Volumes/Expansion/eBird/SMOD_global/GHSL_filt_3cat.tif")
plot(filt)
filt

extent <- c(-8000000, -6000000, 4000000, 6000000)
colors=c("#009E73", "#CC79A7", "#000000")
plot(filt, col=colors, ext=extent)

png(file="new.raster.example.png")
plot(filt, col=colors, ext=extent)
dev.off()




############### Running with different thresholds of modification ############

## threshold of 0.25
GHSL.test <- GHSL
GHSL.test[(dat > 0.25 & GHSL.test %in% c(11,12,13))] <- NA # turn everything that is natural and over 0.5 as NA
GHSL.test[(GHSL.test==10)] <- NA # turn water into NA
plot(GHSL.test)
writeRaster(GHSL.test, "/Volumes/Backup/eBird/SMOD_global/GHSL_filtHighThreshold.tif", overwrite=TRUE)


## threshold of 0.75
GHSL.test <- GHSL
GHSL.test[(dat > 0.75 & GHSL.test %in% c(11,12,13))] <- NA # turn everything that is natural and over 0.5 as NA
GHSL.test[(GHSL.test==10)] <- NA # turn water into NA
GHSL[(GHSL==10)] <- NA
plot(GHSL.test)
writeRaster(GHSL.test, "/Volumes/Backup/eBird/SMOD_global/GHSL_filtLowThreshold.tif", overwrite=TRUE)






