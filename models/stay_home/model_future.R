library(tidyverse)
library(magrittr)
library(feather)
library(rstanarm)
options(mc.cores = 4)

## Read data
model <- readRDS("./model.rds")

county_future <- read_feather("../../data/county_future_stayhome.feather")
#dim(county_future)
#summary(county_future$date[county_future$index_desc == 1])

county_future %<>%
  filter(date <= as.Date("2020-05-05"))
#dim(county_future)

#### #### 
## county_future_fit

county_future_fit <- model %>%
  posterior_predict(county_future, draws = 500)

saveRDS(county_future_fit, "./county_future_fit.rds")
