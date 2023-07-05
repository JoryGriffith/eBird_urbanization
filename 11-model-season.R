########### This is the script for modelling the seasonal differences in bird diversity ###################
library(tidyverse)
library(terra)
library(sp)
library(nlme)
library(ncf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
# Load data
dat <- read.csv("season_model_data.csv")
dat$urban<-as.factor(dat$urban)
dat$BIOME <- as.factor(dat$BIOME)

# make another column with just 3 categories
dat <- dat %>% mutate(urban2=ifelse(urban%in% c(11, 12, 13), 1, ifelse(urban==30, 3, 2)))
dat %>% group_by(urban2) %>% summarise(n=n()) # it worked
dat$urban2 <- as.factor(dat$urban2)
dat$abslat <- abs(dat$lat) # absolute latitude

# Plot relationship 

LDG <- ggplot(dat, aes(x=abs(lat), y=total_SR, color=urban2)) +
  geom_point(alpha=0.1, shape=1) +
  geom_smooth(method="lm") +
  labs(x="Absolute Latitude", y="Species Richness")+
  scale_color_manual(values=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw() +
  facet_wrap(~season)# has the expected relationship!

# Plot data coverage
dat.wint  <- dat %>% filter(season=="winter") # 20,942
plot.wint <-  ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat.wint, aes(x=long, y=lat), color="cornflowerblue", size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude")+
  theme_bw()
ggsave(plot.wint, file="winter.coverage.png", height=8, width=6)

dat.sum  <- dat %>% filter(season=="summer") # 16,116
plot.sum <-  ggplot(data=world)+
  geom_sf() +
  geom_point(data=dat.sum, aes(x=long, y=lat), color="orange", size=0.03) +
  coord_sf(crs = 4326, expand = FALSE) +
  scale_color_viridis_c(na.value = NA, option="B")+
  labs(x="Longitude", y="Latitude")+
  theme_bw()
ggsave(plot.sum, file="summer.coverage.png", height=8, width=6)

dat$season <- factor(dat$season, levels = c("winter", "summer"),
                  labels = c("Winter", "Summer"))

dat$urban2 <- factor(dat$urban2, levels = c("1", "2", "3"),
                     labels = c("Natural", "Suburban", "Urban"))
# model
mod1 <- lm(sqrt(total_SR) ~ abslat * urban2 * season + hemisphere +
                   abslat:hemisphere + BIOME + log(number_checklists), dat)
summary(mod1)
anova(mod1) # triple interaction is significant

# look at model residuals
plot(mod1)
# looking pretty good


# run a gls
gls1 <- gls(sqrt(total_SR) ~ abslat * urban2 * season + hemisphere +
              abslat:hemisphere + BIOME + log(number_checklists), dat)


# look at autocorrelation
datSPDF <- dat
coordinates(datSPDF) <- c("long","lat")
plot(gstat::variogram(residuals(gls1) ~
                        1, data = datSPDF, cutoff = 100))
# looks pretty similar

# try directional
plot(gstat::variogram(residuals(gls1, "normalized") ~
                        1, data = datSPDF, cutoff = 100, alpha = c(0, 45, 90, 135)))
# definitely looks autocorrelated

# some plots
library(jtools)
effect_plot(mod1, abslat)
effect_plot(mod1, hemisphere) # higher richness is N
effect_plot(mod1, BIOME)
effect_plot(mod1, number_checklists)
effect_plot(mod1, season) # higher richness in summer, good

library(interactions)
interact_plot(mod1, abslat, urban2, interval=TRUE) # whoa
interact_plot(mod1, abslat, season, interval=TRUE)
interact_plot(mod1, abslat, hemisphere, interval=TRUE) # interesting, it is steeper in the southern hemisphere
# steeper gradient in winter as predicted

mod.test <- lm(sqrt(total_SR) ~ abslat * urban2 * season, dat)
preds <- (predict(mod.test))^2 # predict and back transform
# plot model predictions
ggplot(dat, aes(x=abslat, y=sqrt(total_SR), color=urban2)) + 
  facet_wrap(~season)+
  geom_line(data = cbind(dat, preds), aes(y = preds)) # it is looking fucked up because there are so many other variables that i am not including in the plot

summary(mod1)
anova(mod1)
#  geom_point(alpha=0.1)
library(emmeans)
library(ggeffects)

emmeans(mod1, specs=c("abslat", "urban2", "season")) 
lstrends(mod1, pairwise ~ urban2, var="abslat", at=c(season="Winter")) # comapre slopes in winter
lstrends(mod1, pairwise ~ urban2, var="abslat", at=c(season="Summer")) # compare slopes in summer

predicted <- ggpredict(mod1, terms = c("abslat", "urban2", "season")) # looks the same whether sqrt included in model or not

seasonal.results.plot<-
  plot(predicted, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#CC79A7", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.text=element_text(size=22), legend.spacing.y = unit(1, 'cm'), legend.position="none")+
  guides(fill = guide_legend(byrow = TRUE))
seasonal.results.plot
ggsave(seasonal.results.plot, file="seasonal.results.plot.png", height=6, width=12)
ggsave(seasonal.results.plot, file="seasonal.results.plot.conclusion.png", height=6, width=7)

predicted2 <- predicted %>% filter(group=="Natural")
seasonal.results.plot2<-
  plot(predicted2, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#009E73")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.text=element_text(size=22), legend.spacing.y = unit(1, 'cm'), legend.position="none")+
  guides(fill = guide_legend(byrow = TRUE))
seasonal.results.plot2
ggsave(seasonal.results.plot2, file="seasonal.results2.plot.png", height=6, width=12)
ggsave(seasonal.results.plot2, file="seasonal.results2.plot.conclusion.png", height=6, width=7) # less wide for conclusion

predicted3 <- predicted %>% filter(group%in%c("Natural", "Urban"))
seasonal.results.plot3<-
  plot(predicted3, facet = TRUE, add.data=TRUE, dot.size=0.5, alpha=0.4, dot.alpha=0.3, 
       line.size=1.5, show.title=FALSE, colors=c("#009E73", "#000000")) +
  theme_bw()+
  labs(x="Absolute Latitude", y="Species Richness", color="Urban")+
  theme(text=element_text(size=20), legend.text=element_text(size=22), legend.spacing.y = unit(1, 'cm'), legend.position="none")
seasonal.results.plot3
ggsave(seasonal.results.plot3, file="seasonal.results.plot3.png", height=6, width=12)
ggsave(seasonal.results.plot3, file="seasonal.results3.plot.conclusion.png", height=6, width=7)
# it back transforms automatically but does not transform the errors, might need to fix later



# subsample and run model to look at correlation
dat.samp <- dat[sample(nrow(dat), 1000), ]

gls1.samp <- gls(sqrt(total_SR) ~ abslat * urban2 * season + hemisphere +
              abslat:hemisphere + BIOME + log(number_checklists), dat.samp)
summary(gls1.samp)
anova(gls1.samp)
# look at autocorrelation
dat.samp$gls1.samp <- residuals(gls1.samp)
hist(dat.samp$gls1.samp)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$gls1.samp, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))

# try to add a correlation structure
gls1.samp.cor <- update(gls1, correlation=corSpatial(form = ~ lat + long | season , nugget=TRUE)) 
beep()
# need to add a grouping factor for season - some overlapping points in different seasons
?gls
dat.samp$gls1.samp.cor <- residuals(gls1.samp.cor)
residsI <- spline.correlog(x=dat.samp$long, y=dat.samp$lat, z=dat.samp$gls1.samp.cor, resamp=100, quiet=TRUE) 
plot(residsI,xlim=c(0,20))
AIC(gls1.samp, gls1.samp.cor)



