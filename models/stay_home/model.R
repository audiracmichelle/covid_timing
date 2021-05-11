library(tidyverse)
library(magrittr)
library(feather)
library(rstanarm)
options(mc.cores = 4)

## Read county_train
county_train <- read_feather("../../data/county_train_stayhome.feather")
#length(unique(county_train$fips))

## Train model
model = stan_glmer.nb(
  y ~
    poly(days_since_thresh, 2) * (nchs + college + age_65_plus + black + hispanic) + 
    (poly(days_since_thresh, 2) | fips) +
    days_since_intrv_stayhome:intrv_stayhome + 
    I(days_since_intrv_stayhome^2):intrv_stayhome + 
    days_since_intrv_stayhome:intrv_stayhome:days_btwn_stayhome_thresh +
    I(days_since_intrv_stayhome^2):intrv_stayhome:days_btwn_stayhome_thresh + 
    days_since_intrv_stayhome:intrv_stayhome:(nchs) +
    I(days_since_intrv_stayhome^2):intrv_stayhome:(nchs)    
  ,
  offset = log(pop),
  data=county_train,
  # algorithm="meanfield",
  iter = 4000,
  warm = 3750,
  chains=4,
  # adapt_iter = 2500,
  QR=TRUE
)

saveRDS(model, paste("./model.rds", sep = ""))

model = readRDS(paste("./model.rds", sep = ""))

#### #### 
## county_fit

county_fit <- model %>%
  posterior_predict(county_train, draws = 500)

county_lp = posterior_linpred(model, newdata=county_train, draws = 500)

saveRDS(county_fit, "./county_fit.rds")
saveRDS(county_lp, "./county_fit_lp.rds")
