library(tidyverse)
library(magrittr)
library(feather)
library(rstanarm)
options(mc.cores=4)

## Read county_train
county_train <- read_feather("../../data/county_train_decrease.feather")
#length(unique(county_train$fips))

## Train model
model = stan_glmer.nb(
  y ~
    poly(days_since_thresh, 2) * (nchs + college + age_65_plus + black + hispanic) + 
    (poly(days_since_thresh, 2) | fips) +
    days_since_intrv_decrease:intrv_decrease + 
    I(days_since_intrv_decrease^2):intrv_decrease + 
    days_since_intrv_decrease:intrv_decrease:days_btwn_decrease_thresh +
    I(days_since_intrv_decrease^2):intrv_decrease:days_btwn_decrease_thresh + 
    days_since_intrv_decrease:intrv_decrease:nchs +
    I(days_since_intrv_decrease^2):intrv_decrease:nchs 
  ,
  offset = log(pop),
  data=county_train,
  # algorithm="meanfield",
  iter = 2500,
  warm = 2250,
  chains = 4,
  # adapt_iter = 2500,
  QR=TRUE
)

# saveRDS(model, paste("./model.rds", sep = ""))
model = readRDS(paste("./model.rds", sep = ""))

#### #### 
## county_fit

county_fit <- model %>%
  posterior_predict(county_train, draws = 500)
county_fit_lp <- model %>%
  posterior_linpred(newdata=county_train, draws = 500)

saveRDS(county_fit, "./county_fit.rds")
saveRDS(county_fit, "./county_fit_lp.rds")
