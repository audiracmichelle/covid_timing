library(tidyverse)
library(magrittr)
library(feather)
library(collections)
library(igraph)
library(lubridate)
options(mc.cores = parallel::detectCores())
source("../utils.R")


## Read county_train
county_train <- read_feather("../data/county_train_.feather") %>%   # 1454 counties
  filter(date <= ymd("20200420")) %>%   # 1021 counties
  group_by(fips) %>%
  filter(
    max(days_since_thresh) >= 7,  # min data points, 909 counties
    max(cum_deaths) >= 1 # there as an outbreak, 400 counties
  ) %>%
  # mutate(mask = !((state == "New York") & (days_since_intrv_decrease >= 12))) %>%
  ungroup()
length(unique(county_train$fips))


model_data = stan_input_data(county_train, type="decrease", lag=0)

model = rstan::stan_model("../stan_models/3b_cable_bent_concave.stan")
print(paste("Compiled model:", "decerase/3b_cable_bent_concave"))

# first run with variational inference,
# it should take ~ 15 mins for default tol (0.01) and ~40min for 50k iters
parnames = c(
  "nchs_pre", "nchs_post", "beta_covars_pre",
  "beta_covars_post", "beta_covars_post",
  "baseline_pre", "baseline_post",
  "overdisp",
  "rand_eff_lin", "state_eff_lin",
  "rand_eff_quad", "state_eff_quad",
  "rand_eff", "state_eff",
  "Omega_rand_eff", "Omega_state_eff",
  "scale_state_eff", "scale_rand_eff",
  "duration_unc", "lag_unc"
)

fit = rstan::vb(
  model, 
  data=model_data,
  adapt_engaged=FALSE,
  eta = 0.25,
  pars=parnames,
  iter=15000,
  tol_rel_obj=0.003,
  adapt_iter=250,
  init="0",
  output_samples=250
)
saveRDS(fit, "decrease_fitted/cable_bent.rds")

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

# now pass solution
fit2 = rstan::sampling(
  model,
  data=model_data,
  chains=nchains,
  iter=5000,
  warmup=4000,
  save_warmup=FALSE,
  pars=parnames,
  thin=10
  init=init_lists
)

# revised_0 uses the joint dataset
# saveRDS(fit, paste("./model_full_rstan_var_revised_0.rds", sep = ""))

# revised_2 uses the joint dataset with min cum deahts >= 1
# saveRDS(fit, paste("./model_full_rstan_var_revised_2.rds", sep = ""))

# 14 experiment
saveRDS(fit2, "decrease_fitted/cable_bent_mcmc.rds")
# fit = read_rds("models/cable_bent.rds")

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

# duration = rstan::extract(fit, pars="duration")$duration
# hist(duration, col=alpha("blue", 0.5), main="duration")
# lag = rstan::extract(fit, pars="lag")$lag
# hist(lag, col=alpha("blue", 0.5), main="lag")
# meandur = mean(duration)


# # let's validate for some location and then call it a day
# # it's working !
# county_lp_var = exp(rstan::extract(fit, pars="log_rate")$log_rate)
# f1 = "06037"  #L.A
# # f1 = "36081"  # queens NY
# f1 = "53033"  # king county WA
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
# dbtwn = dbtwn[dbtwn$days_since_intrv_decrease >= 0, ]
# dbtwn = dbtwn$days_since_thresh[1]
# abline(v=dbtwn + 14, lty=2, col="gray")
# abline(v=dbtwn + 14 + meandur, lty=3, col="gray")
# abline(v=dbtwn + 14 - meandur, lty=3, col="gray")
# title(sprintf("FIPS %s", f1))


# predicted = my_posterior_predict_cable_bent(fit, county_train, type="decrease", eval_pre=FALSE)
# pre_term = apply(predicted$pre_term[ ,ix], 2, median)
# post_term = apply(predicted$post_term[ ,ix], 2, median)
# log_yhat = apply(predicted$log_yhat[, ix], 2, median)

# # predict counterfactual for up

# up = 10
# evaldata_up = county_train %>% 
#   mutate(
#     days_since_intrv_decrease = days_since_intrv_decrease - up,
#     days_btwn_decrease_thresh = days_btwn_decrease_thresh + up,
#   )

# predicted_up = my_posterior_predict_cable_bent(fit, evaldata_up, type="decrease")
# log_yhat_up = apply(predicted_up$log_yhat[, ix], 2, median)

# down = 10
# evaldata_down = county_train %>% 
#   mutate(
#     days_since_intrv_decrease = days_since_intrv_decrease + down,
#     days_btwn_decrease_thresh = days_btwn_decrease_thresh - down,
#   )

# predicted_down = my_posterior_predict_cable_bent(fit, evaldata_down, type="decrease")
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

# soverdisp = rstan::extract(fit, pars="overdisp")$overdisp
# hist(overdisp, col=alpha("blue", 0.5), main="overdisp posterior")

