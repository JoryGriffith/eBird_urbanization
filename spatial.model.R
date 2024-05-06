##### Spatial analysis script ##########
library(spdep)

total.dat <- read.csv("modeling_data.csv")

full.model <- lm(sqrt(total_SR) ~ abslat * urban2 * hemisphere + 
                   precip + log(number_checklists) + elevation, total.dat)

dat.sf <- st_as_sf(dat, coords=c("long", "lat")) 
nb.dist.band <- dnearneigh(dat.sf, 0, 2000) # band of neighbors within 200 km
distances <- nbdists(nb.dist.band,dat.sf) # get distances
invd1 <- lapply(distances, function(x) (1/x)) # get inverse distances
invd.weights <- nb2listw(nb.dist.band, glist = invd1, style = "B") # turn into nblistw



# construct autocovariate
ac <- autocov_dist(resid(full.model), dat.sf, nbs=200, type="inverse")

# add this to model
glm.spat <- glm(total_SR ~ abslat * urban2 + log(number_checklists) + ac, data=dat.samp, family=poisson(link="log"))

AIC(full.model, glm.spat)
# spatial much better
# now check moran's I
moran <- moran.test(residuals.glm(glm.spat), invd.weights)

saveRDS(glm.spat, file="spatial.model.rds")









