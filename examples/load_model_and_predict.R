library(tidyverse)
library(rstan)
library(feather)
options(mc.cores = parallel::detectCores())
source("./utils.R")

## Note:
# Currently, only the fitted stan model is saved without the data
# that was used to train it. This has the disadvantage of having to reload
# the data again when doing predictions (Step 1 below).
# In the future, the data along with the fitted model 
# to save this extra step. 

## 1. Load dataset and generate stan model data

county_train <- read_feather("../data/county_train_.feather") %>%   # 1454 counties
  filter(date <= ymd("20200420")) %>%   # 1021 counties
  group_by(fips) %>%
  filter(
    max(days_since_thresh) >= 7,  # min data points, 909 counties
    max(cum_deaths) >= 1 # there as an outbreak, 400 counties
  ) %>%
  ungroup()
model_data = stan_input_data(county_train, lag=14, type="decrease")

## 2. Load fitted model

fit = read_rds("fit/decrease_fitted/7_states.rds")

## 3. Predict

# The posterior_predict function returns 4 objects:
# - log_yhat  (linear predictor)
# - rand_eff  (random effects component)
# - pre_term  (the pre-intervention curve without the random effects
#              including demographic covariates and nchs interactions)
# - post_term (the post-intevention term (bending of the curve) including
#              effect by timing of intervention and nchs interaction)
# - y_samples (these are simulations of the posterior of the prediction
#              useful for generating confidence intervals that include 
#              the observational noise. If only interested in effects 
#              is is better to use yhat = exp(log_yhat) which removes the observational
#              noise, or yhat_no_re = exp(log_yhat - rand_eff)) which also removes
#              the variation due to random effects.

# Here is how to call the function as well as the defaults

results = posterior_predict(
  fit,
  model_data,
  shift_timing = 0,  # change to positive/negative for delays/anticipation  (resp.) of intervention timing
  states = TRUE,  # include for models that use states random effects
  rand_lag = FALSE,  # if learning lag
  cable_bent = FALSE,  # if learning lag and smooth transition duration
  termporal = FALSE  # if model has an AR(1) component
)