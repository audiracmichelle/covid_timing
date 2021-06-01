library(tidyverse)
library(magrittr)
library(feather)
library(collections)
library(igraph)
library(lubridate)
options(mc.cores = parallel::detectCores())
source("../utils.R")

## exclude_ny?
exclude_ny = FALSE

## Read county_train
county_train <- read_feather("../data/county_train_.feather") %>%   # 1454 counties
  filter(date <= ymd("20200420")) %>%   # 1021 counties
  group_by(fips) %>%
  filter(
    max(days_since_thresh) >= 7,  # min data points, 909 counties
    max(cum_deaths) >= 1 # there as an outbreak, 400 counties
  ) %>%  
  # mutate(mask = !((state == "New York") & (days_since_intrv_decrease >= 14))) %>%
  ungroup()
length(unique(county_train$fips))


model_data = stan_input_data(county_train, type="decrease", lag=0, use_mask=exclude_ny)
model = rstan::stan_model("../stan_models/12b_cable_bent_ar1_concave.stan")
print(paste("Compiled model:", "decrease/12b_cable_bent_ar1_concave"))

parnames = c(
  "nchs_pre", "nchs_post", "beta_covars_pre",
  "beta_covars_post", "beta_covars_post",
  "baseline_pre", "baseline_post",
  "overdisp",
  "rand_eff_lin", "state_eff_lin",
  "rand_eff_quad", "state_eff_quad",
  "state_eff", "rand_eff",
  "Omega_rand_eff", "Omega_state_eff",
  "scale_state_eff", "scale_rand_eff",
  "duration_unc", "lag_unc", "autocor_unc", "time_term"
)

fit = rstan::vb(
  model, 
  data=model_data,
  adapt_engaged=FALSE,
  eta = 0.25,
  iter=50000,
  tol_rel_obj=0.003,
  adapt_iter=250,
  pars=parnames,
  init="0",
  output_samples=250
)


# par
pars = rstan::extract(fit, pars=parnames)

# create list of parameter inialization=
nchains = 4
init_lists = map(1:nchains, function(i) {
  map(pars, function(par) {
    if (length(dim(par))==1)
      return (par[i])
    if (length(dim(par))==2)
      return (par[i, ])
    if (length(dim(par))==3)
      return (par[i, , ])
    print("error")
  })
})
# annoying but must do below becaue R indexing kills a dimension
for (i in 1:nchains) {
  init_lists[[i]]$beta_covars_post = matrix(init_lists[[i]]$beta_covars_post, nrow=1)
  init_lists[[i]]$state_eff_quad = matrix(init_lists[[i]]$state_eff_quad, ncol=1)
  init_lists[[i]]$rand_eff_quad = matrix(init_lists[[i]]$rand_eff_quad, ncol=1)
}

saveRDS(fit, "decrease_fitted/12_cable_bent_ar1.rds")
# now pass solution
fit2 = rstan::sampling(
  model,
  data=model_data,
  chains=nchains,
  iter=2500,
  warmup=2000,
  save_warmup=FALSE,
  pars=parnames,
  thin=10
  init=init_lists
)

saveRDS(fit2, "decrease_fitted/12_cable_bent_ar1_mcmc.rds")

# revised_0 uses the joint dataset
# saveRDS(fit, paste("./model_full_rstan_var_revised_0.rds", sep = ""))

# revised_2 uses the joint dataset with min cum deahts >= 1
# saveRDS(fit, paste("./model_full_rstan_var_revised_2.rds", sep = ""))

# 14 experiment
# if (!exclude_ny) {
#   saveRDS(fit, "models/12_cable_bent_ar1.stan")
#   saveRDS(fit2, "models/12_cable_bent_ar1_mcmc.stan")
#   fit = read_rds("models/12_cable_bent_ar1.stan")
# } else {
#   saveRDS(fit2, "models/12_cable_bent_ar1_mcmc_no_ny.stan")
#   fit = read_rds("models/12_cable_bent_ar1.stan")
# }

# removes the full state of ny
# saveRDS(fit, paste("./model_full_rstan_var_revised_no_ny.rds", sep = ""))

# this one uses the old dataset
# saveRDS(fit, paste("./model_full_rstan_var.rds", sep = ""))

# saveRDS(fit2, paste("./model_full_rstan_mcmc.rds", sep = ""))
# fit = readRDS("./model_full_rstan_var.rds")

# model = readRDS(paste("./model_full_rstan.rds", sep = ""))

#### #### 
## county_fit

# county_fit <- model 
#   posterior_predict(county_train, draws = 500)
# county_fit_var = rstan::extract(fit, pars="y_new")$y_new[-(1:500), ]

# county_fit = rstan::extract(fit2, pars="y_new")$y_new
# county_lp = exp(rstan::extract(fit2, pars="log_rate")$log_rate)

# saveRDS(county_fit_var, "./county_fit_var.rds")
# saveRDS(county_lp_var, "./county_fit_lp_var.rds")

# saveRDS(county_fit, "./county_fit.rds")
# saveRDS(county_lp, "./county_fit_lp.rds")
# fit = fit2
# duration = rstan::extract(fit, pars="duration")$duration
# hist(duration, col=alpha("blue", 0.5), main="duration")
# lag = rstan::extract(fit, pars="lag")$lag
# hist(lag, col=alpha("blue", 0.5), main="lag")
# meandur = mean(duration)


# # let's validate for some location and then call it a day
# # it's working !
# county_lp_var = exp(rstan::extract(fit, pars="log_rate")$log_rate)
# f1 = "06037"  #L.A
# f1 = "36081"  # queens NY
# # f1 = "53033"  # king county WA
# ix = which(county_train$fips == f1)

# yi = county_train$y[ix]
# yhati = apply(county_lp_var[ ,ix], 2, mean)
# yhati_95 = apply(county_lp_var[ ,ix], 2, quantile, .95)
# yhati_05 = apply(county_lp_var[ ,ix], 2, quantile, .05)

# # ymeani = apply(county_fit_var[ ,ix], 2, mean)
# # plot(yi, ylim=c(0, 15))
# plot(yi)
# lines(yhati, col="red")
# lines(yhati_95, col="blue", lty=2)
# lines(yhati_05, col="blue", lty=2)
# dbtwn = county_train[ix, ]
# dbtwn = dbtwn[dbtwn$days_since_intrv_stayhome >= 0, ]
# dbtwn = dbtwn$days_since_thresh[1]
# abline(v=dbtwn + 14, lty=2, col="gray")
# abline(v=dbtwn + 14 + meandur, lty=3, col="gray")
# abline(v=dbtwn + 14 - meandur, lty=3, col="gray")
# title(sprintf("FIPS %s", f1))


# predicted = posterior_predict(fit, model_data, cable_bent=TRUE, temporal=TRUE, states=TRUE)
# pre_term = apply(predicted$pre_term[ ,ix], 2, median)
# post_term = apply(predicted$post_term[ ,ix], 2, median)
# log_yhat = apply(predicted$log_yhat[, ix], 2, median)

# # predict counterfactual for up
# predicted_up = posterior_predict(fit, model_data, rand_lag=TRUE, temporal=TRUE, states=TRUE, shift_timing = 10)
# log_yhat_up = apply(predicted_up$log_yhat[, ix], 2, median)
# predicted_down = posterior_predict(fit, model_data, rand_lag=TRUE, temporal=TRUE, states=TRUE, shift_timing = -10)
# log_yhat_down = apply(predicted_down$log_yhat[, ix], 2, median)


# plotdata = tibble(
#   # no_intervention=pre_term,
#   # intervention_effect=post_term,
#   observed=log_yhat,
#   late10days=log_yhat_up,
#   early10days=log_yhat_down,
#   date=county_train$date[ix],
#   data=log(0.1 + yi)
# ) %>% 
#   pivot_longer(-date)



# ggplot(plotdata) +
#   geom_line(aes(x=date, y=exp(value), color=name), data=filter(plotdata, name != "data")) +
#   geom_point(aes(x=date, y=exp(value)), color="black", data=filter(plotdata, name == "data")) +
#   # geom_vline(aes(xintercept=date[1] + dbtwn - 1), color="black", lty=2) +
#   geom_vline(aes(xintercept=date[1] + 14 + dbtwn - 1), color="black", lty=2) +
#   geom_vline(aes(xintercept=date[1] + 14 + dbtwn + meandur- 1), color="black", lty=3) +
#   geom_vline(aes(xintercept=date[1] + 14 + dbtwn - meandur- 1), color="black", lty=3) +
#   theme_minimal() +
#   labs(
#     title=sprintf("FIPS %s", f1),
#     subtitle="Counterfactual with/without intervention"
#   )

# overdisp = rstan::extract(fit, pars="overdisp")$overdisp
# hist(overdisp, col=alpha("blue", 0.5), main="overdisp posterior")

# autocor = rstan::extract(fit, pars="autocor")$autocor
# hist(autocor, col=alpha("red", 0.5), main="autocor posterior")


# # summary(fit, pars="beta_covars_post")
# # stats
# # parameter                     mean       sd      2.5%       25%
# #   beta_covars_post[1,1] -8.3174984 1.307300 -9.798384 -9.358742
# # beta_covars_post[1,2] -0.9327208 5.407652 -9.198557 -5.696318
# # stats
# # parameter                     50%       75%     97.5%
# # beta_covars_post[1,1] -8.649420 -7.663765 -4.983993
# # beta_covars_post[1,2] -0.850236  3.406478  8.457728

# # traceplot(fit, pars=c("beta_covars_post", "beta_covars_baseline"))