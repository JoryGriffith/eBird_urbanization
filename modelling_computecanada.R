################ Script for running models in compute canada #############

# Load packages
library(tidyverse)
library(nlme)
library(ncf)
library(parallel)

####### Run model 1

# Load data
dat <- read.csv("modeling_data.csv")
dat$BIOME <- as.factor(dat$BIOME)
dat$urban <- as.factor(dat$urban)

dat <- dat %>% mutate(urban2=ifelse(urban==11, 1, ifelse(urban==30, 3, 2)))
dat$urban2 <- as.factor(dat$urban2)
dat$abslat <- abs(dat$lat)

# run model
gls1 <- gls(sqrt(total_SR) ~ abslat * urban2 + hemisphere + abslat:hemisphere + BIOME + log(number_checklists), dat)

# all possible correlation structures
csSpher <- corSpher(form=~lat+long,nugget=TRUE) # sperical
csExp <- corExp(form=~lat+long,nugget=TRUE) # exponential
csGaus <- corGaus(form=~lat+long,nugget=TRUE) # gaussian
csLin <- corLin(form=~lat+long,nugget=TRUE) # linear
csRatio <- corRatio(form=~lat+long,nugget=TRUE) # ratio

# Update models with possible correlation structures
gls1Spher <- update(gls1, correlation=csSpher)
gls1Exp <- update(gls1, correlation=csExp)
gls1Gaus <- update(gls1, correlation=csGaus)
gls1Lin <- update(gls1, correlation=csLin)
gls1Ratio <- update(gls1, correlation=csRatio)

# save as RDS
saveRDS(gls1Spher, "gls1.spher.RDS")
saveRDS(gls1Exp, "gls1.exp.RDS")
saveRDS(gls1Gaus, "gls1.gaus.RDS")
saveRDS(gls1Lin, "gls1.lin.RDS")
saveRDS(gls1Ratio, "gls1.ratio.RDS")












