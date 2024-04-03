#### Testing some generalized linear models

library(tidyverse)
library(ggeffects)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(marginaleffects)
library(nlme)


total.dat <- read.csv("modeling_data.csv")

# linear model with all the variables
full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, total.dat)
plot(ggeffects::ggpredict(full.model, terms = c("abslat", "urban2")), facets = FALSE, dot.alpha=0.2) 
# linear model with less variables
simpl.model <- lm(sqrt(total_SR) ~ abslat * urban2 + log(number_checklists), total.dat)
summary(simpl.model)
# r2 is really not reduced much compared to the full model
plot(ggeffects::ggpredict(simpl.model, terms = c("abslat", "urban2")), facets = FALSE, dot.alpha=0.2) 
# gradients are a bit less steep, but results are honestly very similar


# trying glm with a square root link
glm.mod1 <- glm(total_SR ~ abslat * urban2 + log(number_checklists), data=total.dat, family=poisson(link="log"))
AIC(glm.mod1) # the linear model is a lot better
AIC(simpl.model)
#1398967
plot(ggeffects::ggpredict(glm.mod1, terms = c("abslat", "urban2")), facets = FALSE, dot.alpha=0.2) 
summary(glm.mod1)
plot(glm.mod1)
plot(simpl.model)
hist(resid(glm.mod1))
hist(resid(simpl.model))

# trying glm with all the variables
glm.mod2 <- glm(total_SR ~ abslat * urban2 * hemisphere + precip + log(number_checklists) +
                  elevation, data=total.dat, family=poisson(link="log"))

AIC(full.model) # linear model is still a lot better
AIC(glm.mod2)


res <- glm.mod1$residuals

# extract coordinates corresponding to each residual and combine
res <- data.frame(Residuals = res, x = total.dat$long, y = total.dat$lat)

plot <- ggplot(res, aes(x = x, y = y)) + geom_point(aes(colour = Residuals)) + 
  geom_point(shape = 1, aes(colour = sign) ,colour = "black") + 
  scale_colour_gradient2(low = "#E8F3EB", mid = "#FF1C55",
                         high = "#560000", midpoint = 8, space = "Lab",
                         na.value = "grey50", guide = "colourbar")
plot


dat.samp <- total.dat[sample(nrow(total.dat), 20000), ]

glm.mod.samp <- glm(total_SR ~ abslat * urban2 + log(number_checklists), data=dat.samp, family=poisson(link="log"))

library(sf)
library(spdep)
dat.samp.sf <- st_as_sf(dat.samp, coords=c("long", "lat")) 
nb.dist.band <- dnearneigh(dat.samp.sf, 0, 2000) # band of neighbors
distances <- nbdists(nb.dist.band,dat.samp.sf) # get distances
invd1 <- lapply(distances, function(x) (1/x)) # get inverse distances
invd.weights <- nb2listw(nb.dist.band, glist = invd1, style = "B") # turn into nblistw
class(invd1)
moran.test(residuals.glm(glm.mod.samp), invd.weights)

# construct autocovariate
ac <- autocov_dist(resid(glm.mod.samp), dat.samp.sf, nbs=80, type="inverse")
ac
# add this to model
glm.spat <- glm(total_SR ~ abslat * urban2 + log(number_checklists) + ac, data=dat.samp, family=poisson(link="log"))

AIC(glm.mod.samp, glm.spat)
# spatial much better
# now check moran's I
moran.test(residuals.glm(glm.spat), invd.weights)

plot(ggeffects::ggpredict(glm.mod.samp, terms = c("abslat", "urban2")), facets = FALSE, dot.alpha=0.2) 
plot(ggeffects::ggpredict(glm.spat, terms = c("abslat", "urban2")), facets = FALSE, dot.alpha=0.2) 


# trying a linear mixed effects model
library(nlme)
lme1 <- lme(sqrt(total_SR) ~ abslat * urban2 * hemisphere, 
            random = ~ 1 | log(number_checklists) + 1 | precip + 1| elevation, data=total.dat, na.action=na.omit)
AIC(lme1) # mixed effects is worse
AIC(full.model)

summary(lme1)
plot(lme1)



