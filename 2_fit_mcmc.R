library(tidyverse)
library(magrittr)
library(feather)
library(igraph)
library(lubridate)
library(argparse)
library(yaml)
library(reticulate)
library(rstan)
library(abind)

options(mc.cores = 4)
source("./utils.R")

parser = ArgumentParser()
parser$add_argument("--intervention", type="character",
    default="decrease",
    choices=c("decrease", "stayhome"),
    help="intervention_type"
)
parser$add_argument("--bent_cable", action="store_true", default=FALSE,
    help="Fits the model with random lag")
parser$add_argument("--exclude_ny", action="store_true", default=FALSE,
    help="Exclude the counties from NY city")
parser$add_argument("--exclude_cities_post_inter", action="store_true", default=FALSE,
    help="Exclude the intervention data for Wayne, NYC and LA")
parser$add_argument("--lag", type="double", default=14.0, 
    help="Lag to use (time from infection to death).")
parser$add_argument("--dir", type="character", default="", 
    help="Folder where to save output model.")
parser$add_argument("--no_spatial", action="store_false", default=TRUE, dest="spatial",
    help="Adds an intrinsic spatial ICAR term")
parser$add_argument("--no_temporal", action="store_false", default=TRUE, dest="temporal",
    help="Adds am intrinsic tempoeral autoregressive term. The autocorrelation is set by --autocor")
parser$add_argument("--autocor", type="double", default=0.8, 
    help="Autocorrelation for intrinsic term. If autocor=1.0 it uses a random walk.")
parser$add_argument("--spatial_scale", type="double", default=0.0, dest="spatial_scale_fixed",
    help="When positive it fixes the scale of the spatial process. If 0.0 the scale is learned with BYM scaling.")
parser$add_argument("--ar_scale", type="double", default=0.0, dest="ar_scale_fixed",
    help="When positive it fixes the scale of the temporal process. If 0.0 the scale is learned with AUTO SCALING.")
parser$add_argument("--duration", type="double", default=0.0, dest="duration_fixed",
    help="Only used for bent cable model. If different from 0 duration is fixed.")
parser$add_argument("--use_post_inter", action="store_true", default=FALSE,
    help="Use interaction variables for post-trend")
parser$add_argument("--no_pre_inter", action="store_false", default=TRUE, dest="use_pre_inter",
    help="Use interaction variables for pre-trend")
parser$add_argument("--pre_vars", type="character",
    default="college age_65_plus black hispanic",
    help="Control variables for the pre-trend")
parser$add_argument("--pre_inter_vars", type="character",
    default="college age_65_plus black hispanic",
    help="Control variables for the pre-trend  that interact with NCHS")
parser$add_argument("--post_vars", type="character",
    default="college age_65_plus black hispanic",
    help="Control variables for the post-trend")
parser$add_argument("--post_inter_vars", type="character",
    default="",
    help="Control variables for the post-trend that interact with NCHS. The timing (days between intervention and threshold) is always added to the list.")
parser$add_argument("--ar_tight_prior_scale", action="store_true", default=FALSE,
    help="Deprecated!")
parser$add_argument("--nchains", type="integer", default=3, 
    help="Number of chains")
parser$add_argument("--thin", type="integer", default=10, 
    help="Thinning. Adjust to with chains and iters to more less end up with 800~1000 samples. (More will break memory in analysis.)")
parser$add_argument("--iter", type="integer", default=4000, 
    help="Number of iterations for the chain")
parser$add_argument("--warmup", type="integer", default=2000, 
    help="Weramup iterations for MCMC Use a larger number for cold starts.")
parser$add_argument("--init", type="character", default="random",
    help="If distinct from random or '0', it must be the the path to a .rds file with an instance of a fitted model to sample the initial distribution of parameters for the chain.")


args = parser$parse_args()
for (i in 1:length(args)) {
    varname = names(args)[i]
    value = args[[i]]
    print(sprintf("%s: %s", varname, value))
    if (varname %in% c("pre_vars", "post_vars", "pre_inter_vars", "post_inter_vars"))
        value = strsplit(value, " ")[[1]]
    if (varname == "dir")
        value = paste0("results/", value)
    assign(varname, value, envir = .GlobalEnv)
}
dir.create(dir, recursive=TRUE)
write_yaml(args, sprintf("%s/config.yaml", dir))


ny_counties = c(
  "36081",  # queens
  "36085",  # richmond
  "36061",  # new york
  "36047",  # kings
  "36005"  # bronx
)

out_cities = c(
    "Queens (NYC)"="36081",
    "Wayne (MI)"="26163",
    "Los Angeles (CA)"="06037", 
    "Kings (WA)"="53033"
)

## Read c

## Read county_train
county_train <- read_feather("data/county_train_.feather") %>%   # 1454 counties
  group_by(fips) %>%
  filter(date <= ymd("20200420")) %>%   # 1021 counties
  filter(
    max(days_since_thresh) >= 7,  # min data points, 909 counties
    max(cum_deaths) >= 1 # there as an outbreak, 400 counties
  ) %>%
  mutate(interv_type = .env$intervention) %>%
  mutate(days_since_intrv = if_else(interv_type == "decrease", days_since_intrv_decrease, days_since_intrv_stayhome)) %>%
  ungroup()
valid_fips = unique(county_train$fips)


if (exclude_ny)
    county_train$mask = !(county_train$fips %in% ny_counties)
if (exclude_cities_post_inter)
    county_train$mask = !((county_train$fips %in% out_cities) & (county_train$days_since_intrv > lag))

use_mask = (exclude_ny || exclude_cities_post_inter)
edges = read_csv("data/edges_with_knn2_filt100.csv")

model_data = stan_input_data(
  county_train,
  type=intervention,
  lag=lag,
  use_mask=use_mask,
  pre_vars=pre_vars,
  post_vars=post_vars,
  post_inter_vars=post_inter_vars,
  pre_inter_vars=pre_inter_vars,
  use_pre_inter=use_pre_inter,
  use_post_inter=use_post_inter,
  autocor=autocor,
  edges=edges,
  ar_tight_prior_scale=ar_tight_prior_scale,
  spatial_scale_fixed=spatial_scale_fixed,
  bent_cable=bent_cable,
  spatial=spatial,
  temporal=temporal,
  duration=duration_fixed,
  ar_scale_fixed=ar_scale_fixed
)
saveRDS(model_data, paste0(dir, "/model_data.rds"))

pars = c(
  "nchs_pre",
  "nchs_post",
  "beta_covars_pre",
  "beta_covars_pre_inter",
  "beta_covars_post",
  "beta_covars_post_inter",
  "baseline_pre",
  "baseline_post",
  "overdisp",
  "rand_eff_lin",
  "rand_eff_quad",
  "rand_eff",
  "state_eff_lin",
  "state_eff_quad",
  "state_eff",
  "spatial_eff_lin",
  "spatial_eff_quad",
  "spatial_eff",
  "time_term",
  "Omega_rand_eff",
  "Omega_state_eff",
  "scale_state_eff",
  "scale_rand_eff",
  "spatial_scale",
  "ar_scale"
#   "log_rate"
)
if (bent_cable)
    pars = c(pars, c("lag_unc", "duration_unc"))


if (init != "random" && init != "0") {
  prevfit = read_rds(init)
  prevpars = rstan::extract(prevfit, pars=pars)
  init = map(1:nchains, function(i) {
    map(prevpars, function(par) {
      if (length(dim(par))==1)
        return (par[i])
      if (length(dim(par))==2)
        return (par[i, ])
      if (length(dim(par))==3)
        return (par[i, , ])
      print("error")
    })
  })

  # annoying but must do below bacause R indexing kills a dimension (drop=TRUE not helping) when D_post has dimension 1
  for (i in 1:nchains) {
    if (model_data$D_post == 1)
      init[[i]]$beta_covars_post = matrix(init[[i]]$beta_covars_post, nrow=1)
    if (model_data$D_post_inter == 1)
      init[[i]]$beta_covars_post_inter = matrix(init[[i]]$beta_covars_post_inter, ncol=1)
    init[[i]]$rand_eff_quad = matrix(init[[i]]$rand_eff_quad, ncol=1)
    init[[i]]$state_eff_quad = matrix(init[[i]]$state_eff_quad, ncol=1)
    init[[i]]$spatial_eff_quad = matrix(init[[i]]$spatial_eff_quad, ncol=1)
  }
}

if (!bent_cable) {
    model = rstan::stan_model("stan/1_spatiotemporal.stan")
} else{
    model = rstan::stan_model("stan/1_spatiotemporal_bentcable.stan")
}

fit_mcmc = rstan::sampling(
  model,
  data=model_data,
  chains=nchains,
  iter=iter,
  warmup=warmup,
  init=init,
  pars=c("log_rate", pars),
  thin=thin,
  refresh=10 # to measure sampling speed
)
saveRDS(fit_mcmc, paste0(dir, "/fit.rds"))
