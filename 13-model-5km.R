### Modelling with the scaled up data ######

library(tidyverse)
library(terra)

## Load summary data
summary <- read.csv("global_richness_summary_5km.csv")

### Use coverage to threshold the data
coverage <- read.csv("coverage_top500_5km.csv")

quantile(coverage$sampsize_95, 0.95) # 102

# Remove data with a sample size less than 102
summary_filt95 <- summary %>% filter(number_checklists >= 102) # went from 569861 to 57873

#### Now I need to load the urbanization data and figure out a threshold for that

