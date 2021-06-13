library(tidyverse)
library(magrittr)
library(feather)
library(igraph)
library(lubridate)
library(argparse)
library(yaml)
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
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
parser$add_argument("--exclude_cities_post", action="store_true", default=FALSE,
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
parser$add_argument("--iter", type="double", default=100000, 
    help="Max number of iterations for stan variational algorithm")
parser$add_argument("--samples", type="integer", default=250, 
    help="Number of sample samples for rstan vb algorithm")
parser$add_argument("--rel_tol", type="double", default=0.001, 
    help="Relative tolerance for rstan vb algorithm")
parser$add_argument("--eta", type="double", default=0.2, 
    help="eta for variational step (see rstan::vb documentation)")


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


ny_counties = c(
  "36081",  # queens
  "36085",  # richmond
  "36061",  # new york
  "36047",  # kings
  "36005"  # bronx
)

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
  mutate(days_btwn = if_else(interv_type == "decrease", days_btwn_decrease_thresh, days_btwn_stayhome_thresh)) %>%
  ungroup()
valid_fips = unique(county_train$fips)

if (exclude_ny) {
    if (exclude_post_only) {
        county_train$mask = !((county_train$fips %in% ny_counties) & (county_train$days_since_intrv > lag))
    } else {
        county_train$mask = !(county_train$fips %in% ny_counties)
    }
}

edges = read_csv("data/edges_with_knn2_filt100.csv")

model_data = stan_input_data(
  county_train,
  type=intervention,
  lag=lag,
  use_mask=exclude_ny,
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
  ar_scale_fixed=ar_scale_fixed
)

dir.create(dir, recursive=TRUE)
write_yaml(args, sprintf("%s/config.yaml", dir))
file.copy(c("stan/1_spatiotemporal.stan", "utils.R"), dir)

saveRDS(model_data, paste0(dir, "/model_data.rds"))

print("X_pre")
print(head(model_data$X_pre))
print("X_post")
print(head(model_data$X_post))
print("use_pre_inter")
print(model_data$use_pre_inter)
print("use_post_inter")
print(model_data$use_post_inter)
print("autocor")
print(head(model_data$autocor))

if (!bent_cable) {
    model = rstan::stan_model("stan/1_spatiotemporal.stan")
} else{
    model = rstan::stan_model("stan/1_spatiotemporal_bentcable.stan")
}
    

# parameters form the model to save samples from
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
  "Omega_rand_eff",
  "Omega_state_eff",
  "scale_state_eff",
  "scale_rand_eff",
  "spatial_eff_lin",
  "spatial_eff_quad",
  "spatial_eff",
  "time_term",
  "spatial_scale",
  "ar_scale",
  "log_rate"
)

if (bent_cable)
    pars = c(pars, c("lag_unc", "duration_unc"))

# pre-fitting the variational model should take ~ 1h to convergence
fit_vb = rstan::vb(
  model, 
  data=model_data,
  adapt_engaged=FALSE,
  eta = eta,
  iter=iter,
  tol_rel_obj=rel_tol,
  init="0",
  pars=pars,
  output_samples=samples
)
saveRDS(fit_vb, paste0(dir, "/fit.rds"))
