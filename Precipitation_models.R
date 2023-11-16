######## Trying model with precipitation instead of biome #######

library(terra)
library(tidyverse)
library(marginaleffects)
library(patchwork)
library(ggeffects)

#### Load new precipitation data
#precip.test <- rast("/Users/jorygriffith/Downloads/3B-MO-GIS.MS.MRG.3IMERG.20210101-S000000-E235959.01.V07A/3B-MO-GIS.MS.MRG.3IMERG.20210101-S000000-E235959.01.V07A.total.accum.tif")
#plot(precip.test)
#summary(values(precip.test))
## no NAs so that's good, but there are some crazy values
#
## Going to load each month and then average them and then add up to get year accumulation
## Jan
#jan17 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jan/3B-MO-GIS.MS.MRG.3IMERG.20170101-S000000-E235959.01.V07A.total.accum.tif")
#jan18 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jan/3B-MO-GIS.MS.MRG.3IMERG.20180101-S000000-E235959.01.V07A.total.accum.tif")
#jan19 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jan/3B-MO-GIS.MS.MRG.3IMERG.20190101-S000000-E235959.01.V07A.total.accum.tif")
#jan20 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jan/3B-MO-GIS.MS.MRG.3IMERG.20200101-S000000-E235959.01.V07A.total.accum.tif")
#jan21 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jan/3B-MO-GIS.MS.MRG.3IMERG.20210101-S000000-E235959.01.V07A.total.accum.tif")
#jan22 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jan/3B-MO-GIS.MS.MRG.3IMERG.20220101-S000000-E235959.01.V07A.total.accum.tif")
#jan_stack <- c(jan17, jan18, jan19, jan20, jan21, jan22)
#jan_avg <- app(jan_stack, fun=mean)
#plot(jan_avg)
#
## Feb 
#feb17 <- rast("/Users/jorygriffith/Desktop/Precipitation/Feb/3B-MO-GIS.MS.MRG.3IMERG.20170201-S000000-E235959.02.V07A.total.accum.tif")
#feb18 <- rast("/Users/jorygriffith/Desktop/Precipitation/Feb/3B-MO-GIS.MS.MRG.3IMERG.20180201-S000000-E235959.02.V07A.total.accum.tif")
#feb19 <- rast("/Users/jorygriffith/Desktop/Precipitation/Feb/3B-MO-GIS.MS.MRG.3IMERG.20190201-S000000-E235959.02.V07A.total.accum.tif")
#feb20 <- rast("/Users/jorygriffith/Desktop/Precipitation/Feb/3B-MO-GIS.MS.MRG.3IMERG.20200201-S000000-E235959.02.V07A.total.accum.tif")
#feb21 <- rast("/Users/jorygriffith/Desktop/Precipitation/Feb/3B-MO-GIS.MS.MRG.3IMERG.20210201-S000000-E235959.02.V07A.total.accum.tif")
#feb22 <- rast("/Users/jorygriffith/Desktop/Precipitation/Feb/3B-MO-GIS.MS.MRG.3IMERG.20220201-S000000-E235959.02.V07A.total.accum.tif")
#feb_stack <- c(feb17, feb18, feb19, feb20, feb21, feb22)
#feb_avg <- app(feb_stack, fun=mean)
#
## Mar
#mar17 <- rast("/Users/jorygriffith/Desktop/Precipitation/Mar/3B-MO-GIS.MS.MRG.3IMERG.20170301-S000000-E235959.03.V07A.total.accum.tif")
#mar18 <- rast("/Users/jorygriffith/Desktop/Precipitation/Mar/3B-MO-GIS.MS.MRG.3IMERG.20180301-S000000-E235959.03.V07A.total.accum.tif")
#mar19 <- rast("/Users/jorygriffith/Desktop/Precipitation/Mar/3B-MO-GIS.MS.MRG.3IMERG.20190301-S000000-E235959.03.V07A.total.accum.tif")
#mar20 <- rast("/Users/jorygriffith/Desktop/Precipitation/Mar/3B-MO-GIS.MS.MRG.3IMERG.20200301-S000000-E235959.03.V07A.total.accum.tif")
#mar21 <- rast("/Users/jorygriffith/Desktop/Precipitation/Mar/3B-MO-GIS.MS.MRG.3IMERG.20210301-S000000-E235959.03.V07A.total.accum.tif")
#mar22 <- rast("/Users/jorygriffith/Desktop/Precipitation/Mar/3B-MO-GIS.MS.MRG.3IMERG.20220301-S000000-E235959.03.V07A.total.accum.tif")
#mar_stack <- c(mar17, mar18, mar19, mar20, mar21, mar22)
#mar_avg <- app(mar_stack, fun=mean)
#
## Apr
#apr17 <- rast("/Users/jorygriffith/Desktop/Precipitation/Apr/3B-MO-GIS.MS.MRG.3IMERG.20170401-S000000-E235959.04.V07A.total.accum.tif")
#apr18 <- rast("/Users/jorygriffith/Desktop/Precipitation/Apr/3B-MO-GIS.MS.MRG.3IMERG.20180401-S000000-E235959.04.V07A.total.accum.tif")
#apr19 <- rast("/Users/jorygriffith/Desktop/Precipitation/Apr/3B-MO-GIS.MS.MRG.3IMERG.20190401-S000000-E235959.04.V07A.total.accum.tif")
#apr20 <- rast("/Users/jorygriffith/Desktop/Precipitation/Apr/3B-MO-GIS.MS.MRG.3IMERG.20200401-S000000-E235959.04.V07A.total.accum.tif")
#apr21 <- rast("/Users/jorygriffith/Desktop/Precipitation/Apr/3B-MO-GIS.MS.MRG.3IMERG.20210401-S000000-E235959.04.V07A.total.accum.tif")
#apr22 <- rast("/Users/jorygriffith/Desktop/Precipitation/Apr/3B-MO-GIS.MS.MRG.3IMERG.20220401-S000000-E235959.04.V07A.total.accum.tif")
#apr_stack <- c(apr17, apr18, apr19, apr20, apr21, apr22)
#apr_avg <- app(apr_stack, fun=mean)
#
## May 
#may17 <- rast("/Users/jorygriffith/Desktop/Precipitation/May/3B-MO-GIS.MS.MRG.3IMERG.20170501-S000000-E235959.05.V07A.total.accum.tif")
#may18 <- rast("/Users/jorygriffith/Desktop/Precipitation/May/3B-MO-GIS.MS.MRG.3IMERG.20180501-S000000-E235959.05.V07A.total.accum.tif")
#may19 <- rast("/Users/jorygriffith/Desktop/Precipitation/May/3B-MO-GIS.MS.MRG.3IMERG.20190501-S000000-E235959.05.V07A.total.accum.tif")
#may20 <- rast("/Users/jorygriffith/Desktop/Precipitation/May/3B-MO-GIS.MS.MRG.3IMERG.20200501-S000000-E235959.05.V07A.total.accum.tif")
#may21 <- rast("/Users/jorygriffith/Desktop/Precipitation/May/3B-MO-GIS.MS.MRG.3IMERG.20210501-S000000-E235959.05.V07A.total.accum.tif")
#may22 <- rast("/Users/jorygriffith/Desktop/Precipitation/May/3B-MO-GIS.MS.MRG.3IMERG.20220501-S000000-E235959.05.V07A.total.accum.tif")
#may_stack <- c(may17, may18, may19, may20, may21, may22)
#may_avg <- app(may_stack, fun=mean)
#
## Jun3
#jun17 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jun/3B-MO-GIS.MS.MRG.3IMERG.20170601-S000000-E235959.06.V07A.total.accum.tif")
#jun18 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jun/3B-MO-GIS.MS.MRG.3IMERG.20180601-S000000-E235959.06.V07A.total.accum.tif")
#jun19 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jun/3B-MO-GIS.MS.MRG.3IMERG.20190601-S000000-E235959.06.V07A.total.accum.tif")
#jun20 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jun/3B-MO-GIS.MS.MRG.3IMERG.20200601-S000000-E235959.06.V07A.total.accum.tif")
#jun21 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jun/3B-MO-GIS.MS.MRG.3IMERG.20210601-S000000-E235959.06.V07A.total.accum.tif")
#jun22 <- rast("/Users/jorygriffith/Desktop/Precipitation/Jun/3B-MO-GIS.MS.MRG.3IMERG.20220601-S000000-E235959.06.V07A.total.accum.tif")
#jun_stack <- c(jun17, jun18, jun19, jun20, jun21, jun22)
#jun_avg <- app(jun_stack, fun=mean)
#
## July
#list_csv_files <- list.files(path = "/Users/jorygriffith/Desktop/Precipitation/Jun", pattern="*.tif")
#
#names <- c("jul17", "jul18", "jul19", "jul20", "jul21", "jul22")
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         rast(paste("/Users/jorygriffith/Desktop/Precipitation/Jun/", list_csv_files[i], sep="")))
#}
#jul_avg <- app(c(jul17, jul18, jul19, jul20, jul21, jul22), fun=mean)
#
#
## August
#list_csv_files <- list.files(path = "/Users/jorygriffith/Desktop/Precipitation/Aug", pattern="*.tif")
#
#names <- c("aug17", "aug18", "aug19", "aug20", "aug21", "aug22")
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         rast(paste("/Users/jorygriffith/Desktop/Precipitation/Aug/", list_csv_files[i], sep="")))
#}
#aug_avg <- app(c(aug17, aug18, aug19, aug20, aug21, aug22), fun=mean)
#
#
## September
#list_csv_files <- list.files(path = "/Users/jorygriffith/Desktop/Precipitation/Sept", pattern="*.tif")
#
#names <- c("sep17", "sep18", "sep19", "sep20", "sep21", "sep22")
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         rast(paste("/Users/jorygriffith/Desktop/Precipitation/Sept/", list_csv_files[i], sep="")))
#}
#sep_avg <- app(c(sep17, sep18, sep19, sep20, sep21, sep22), fun=mean)
#
#
## October
#list_csv_files <- list.files(path = "/Users/jorygriffith/Desktop/Precipitation/Oct", pattern="*.tif")
#
#names <- c("oct17", "oct18", "oct19", "oct20", "oct21", "oct22")
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         rast(paste("/Users/jorygriffith/Desktop/Precipitation/Oct/", list_csv_files[i], sep="")))
#}
#oct_avg <- app(c(oct17, oct18, oct19, oct20, oct21, oct22), fun=mean)
#
#
## November
#list_csv_files <- list.files(path = "/Users/jorygriffith/Desktop/Precipitation/Nov", pattern="*.tif")
#
#names <- c("nov17", "nov18", "nov19", "nov20", "nov21", "nov22")
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         rast(paste("/Users/jorygriffith/Desktop/Precipitation/Nov/", list_csv_files[i], sep="")))
#}
#nov_avg <- app(c(nov17, nov18, nov19, nov20, nov21, nov22), fun=mean)
#
## December
#list_csv_files <- list.files(path = "/Users/jorygriffith/Desktop/Precipitation/Dec", pattern="*.tif")
#
#names <- c("dec17", "dec18", "dec19", "dec20", "dec21", "dec22")
#
#for(i in 1:length(list_csv_files)) {                              # Head of for-loop
#  assign(names[i],                                   # Read and store data frames
#         rast(paste("/Users/jorygriffith/Desktop/Precipitation/Dec/", list_csv_files[i], sep="")))
#}
#dec_avg <- app(c(dec17, dec18, dec19, dec20, dec21, dec22), fun=mean)
#
#
##### put all together and get an annual accumulation!
#annual_accum <- app(c(jan_avg, feb_avg, mar_avg, apr_avg, may_avg, jun_avg, jul_avg, aug_avg, sep_avg, oct_avg, nov_avg, dec_avg), fun=sum)
#
#plot(annual_accum)
#
#

# Load data
dat <- read.csv("modeling_data.csv")

# extract for data points
dat$precip <- as.data.frame(terra::extract(annual_accum, dat[,c(27:28)]))$sum
range(dat$precip)

# Run model with precipitation instead of biome
mod.precip <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                          precip + log(number_checklists) + elevation, dat) # decent number of points with no precipitation values, would need to fix

mod.biome <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
            BIOME + log(number_checklists) + elevation, dat)

AIC(mod.precip, mod.biome) # the model for precipitation has a higher AIC than biome


square <- function(x){
  x^2
} 
predicted.full<-avg_predictions(mod.precip, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))
# want to make this over a more limited set of predictors and somehow make not wiggly
#marginal.full<-plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) # need to figure out how to back transform this

mainLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mainLDGplot





predicted.full<-avg_predictions(mod.biome, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))
# want to make this over a more limited set of predictors and somehow make not wiggly
#marginal.full<-plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) # need to figure out how to back transform this

mainLDGplot2 <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mainLDGplot2


mainLDGplot | mainLDGplot2
# look pretty much the same, slope is  little bit steeper in biome plot

# look at model results
hypothesis_test(mod.precip, c("abslat", "urban2"), test = NULL)

hypothesis_test(mod.biome, c("abslat", "urban2"), test = NULL) # confidence intervals are slightly different but slopes are not


##########################
## Now for seasonal
seasonal.dat <- read.csv("season_modeling_data.csv")

seasonal.dat$precip <- as.data.frame(terra::extract(annual_accum, seasonal.dat[,c(15:16)]))$sum

# run model
mod.precip.season <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                          precip + log(number_checklists) + elevation, seasonal.dat)

mod.biome.season <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                         BIOME + log(number_checklists) + elevation, seasonal.dat)
AIC(mod.precip.season, mod.biome.season) # precip is better (but again need to figure out these NAs)

# Plot results
predicted.season <- avg_predictions(mod.precip.season, by=c("abslat", "urban2", "season"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                       season = c("Summer", "Winter")))

#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(seasonal.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot


predicted.season <- avg_predictions(mod.biome.season, by=c("abslat", "urban2", "season"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                       season = c("Summer", "Winter")))

#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot2 <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(seasonal.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())

seasonLDGplot | seasonLDGplot2
# very similar


hypothesis_test(mod.precip.season, c("abslat", "urban2", "season"), test = NULL)

hypothesis_test(mod.biome.season, c("abslat", "urban2", "season"), test = NULL) 
# same overall results which is good



summary(dat) # there are 3058 NAs for precipitation
summary(seasonal.dat) # 2147 NAs for precipitation

?terra::extract





######## CODE FOR USING WORLDCLIM
# January
#jan <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_01.tif")
#summary(values(jan)) # there are a bunch of missing values in the raster
#plot(jan)
#feb <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_02.tif")
#mar <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_03.tif")
#apr <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_04.tif")
#may <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_05.tif")
#jun <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_06.tif")
#jul <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_07.tif")
#aug <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_08.tif")
#sep <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_09.tif")
#oct <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_10.tif")
#nov <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_11.tif")
#dec <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_prec/wc2.1_30s_prec_12.tif")

#year_stack <- c(jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec)

#mean_prec <- app(year_stack, fun=mean)
#plot(mean_prec)
# looks good






# extract data from raster
#dat$precip <- as.data.frame(terra::extract(mean_prec, dat[,c(27:28)]))$mean






########### Extract data from worldclim bioclimatic variables using dismo package
# Load bio12 (annual precipitation)
#precip <- rast("/Users/jorygriffith/Desktop/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
precip <- rast("/Users/jorygriffith/Desktop/wc2.1_5m_bio/wc2.1_5m_bio_12.tif")

# going to try with the larger resolution
plot(precip, ext=c(-100, -70, 30, 60))
precip
crs(precip)


dat <- read.csv("modeling_data.csv")

# extract for data points
dat$precip <- as.data.frame(terra::extract(precip, dat[,c(27:28)], method="bilinear"))$wc2.1_5m_bio_12
# doing bilinear because some values need to be interpolated

mod.precip <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, dat) # decent number of points with no precipitation values, would need to fix

mod.biome <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                  BIOME + log(number_checklists) + elevation, dat)

AIC(mod.precip, mod.biome) # the model for precipitation has a higher AIC than biome


square <- function(x){
  x^2
} 
predicted.full<-avg_predictions(mod.precip, by=c("abslat", "urban2"), transform=square, 
                                newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban")))
# want to make this over a more limited set of predictors and somehow make not wiggly
#marginal.full<-plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) # need to figure out how to back transform this

mainLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.full, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.full, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
mainLDGplot
# precipitation changed the sloe ever so slightly but nothing crazy



#### Run for season
seasonal.dat <- read.csv("season_modeling_data.csv")

seasonal.dat$precip <- as.data.frame(terra::extract(precip, seasonal.dat[,c(15:16)], method="bilinear"))$wc2.1_5m_bio_12
# only two nan values for seasonal data, great!

# run model
mod.precip.season <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                          precip + log(number_checklists) + elevation, seasonal.dat)

mod.biome.season <- lm(sqrt(total_SR) ~ abslat * urban2 * season * hemisphere + 
                         BIOME + log(number_checklists) + elevation, seasonal.dat)
AIC(mod.precip.season, mod.biome.season) # biome is slightly better

# Plot results
predicted.season <- avg_predictions(mod.biome.season, by=c("abslat", "urban2", "season"), transform=square, 
                                    newdata = datagrid(abslat = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70), urban2=c("Natural", "Suburban", "Urban"),
                                                       season = c("Summer", "Winter")))

#predicted.season<-avg_predictions(season.model, by=c("abslat", "urban2", "season"), transform=square, newdata="mean")
seasonLDGplot <- #plot_predictions(full.model, condition=c("abslat", "urban2"), transform=square, points=0.01) + 
  ggplot()+
  geom_point(seasonal.dat, mapping=aes(x=abslat, y=total_SR, color=urban2), size=0.25, alpha=0.1)+
  geom_line(predicted.season, mapping=aes(x=abslat, y=estimate, color=urban2), lwd=1.5)+
  geom_ribbon(predicted.season, mapping=aes(x=abslat, ymax=conf.high, ymin=conf.low, group=urban2), alpha=0.3)+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  scale_fill_manual(values=c("#009E73", "#CC79A7", "#000000"))+
  labs(x="Absolute latitude", y="Species richness")+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0))+
  facet_wrap(~season)+
  theme(legend.title=element_blank(), legend.position = c(.8, .85), text=element_text(size=15), axis.title=element_blank())
seasonLDGplot
# looks pretty much the same!







