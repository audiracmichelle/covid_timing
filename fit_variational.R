library(tidyverse)
library(magrittr)
library(feather)
library(igraph)
library(lubridate)
library(argparse)
library(yaml)
options(mc.cores = 4)
source("./utils.R")

parser = ArgumentParser()
parser$add_argument("--intervention", type="character",
    default="decrease",
    choices=c("decrease", "stayhome"),
    help="intervention_type"
)
parser$add_argument("--exclude_ny", action="store_true", default=FALSE,
    help="Exclude the five counties from NY city")
parser$add_argument("--exclude_post_only", action="store_true", default=FALSE,
    help="Exclude the five counties from NY only after intervention time")
parser$add_argument("--lag", type="double", default=14.0, 
    help="Lag to use (time from infection to death).")
parser$add_argument("--dir", type="character", default="", 
    help="Folder where to save output model.")
parser$add_argument("--no_spatial", action="store_false", default=TRUE, dest="spatial",
    help="Adds an intrinsic spatial ICAR term")
parser$add_argument("--no_temporal", action="store_false", default=TRUE, dest="temporal",
    help="Adds am intrinsic tempoeral autoregressive term. The autocorrelation is set by --autocor")
parser$add_argument("--autocor", type="double", default=0.7, 
    help="Autocorrelation for intrinsic term. If autocor=1.0 it uses a random walk.")
parser$add_argument("--spatial_scale", type="double", default=0.0, dest="spatial_scale_fixed",
    help="When positive it fixes the scale of the spatial process. If 0.0 the scale is learned with BYM scaling.")
parser$add_argument("--use_post_inter", action="store_true", default=FALSE,
    help="Use interaction variables for post-trend")
parser$add_argument("--no_pre_inter", action="store_false", default=TRUE, dest="use_pre_inter",
    help="Use interaction variables for pre-trend")
parser$add_argument("--pre_vars", type="character",
    default="college age_65_plus black hispanic popdensity",
    help="Control variables for the pre-trend")
parser$add_argument("--pre_inter_vars", type="character",
    default="college age_65_plus black hispanic popdensity",
    help="Control variables for the pre-trend  that interact with NCHS")
parser$add_argument("--post_vars", type="character",
    default="college age_65_plus black hispanic popdensity",
    help="Control variables for the post-trend")
parser$add_argument("--post_inter_vars", type="character",
    default="",
    help="Control variables for the post-trend that interact with NCHS. The timing (days between intervention and threshold) is always added to the list.")
parser$add_argument("--ar_tight_prior_scale", action="store_true", default=FALSE,
    help="When selected the prior scaling is 1/N instead of 1/sqrt(N).")
parser$add_argument("--iter", type="double", default=100000, 
    help="Max number of iterations for stan variational algorithm")
parser$add_argument("--samples", type="integer", default=250, 
    help="Number of sample samples for rstan vb algorithm")
parser$add_argument("--rel_tol", type="double", default=0.001, 
    help="Relative tolerance for rstan vb algorithm")

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

if (exclude_ny) {
    if (exclude_post_only) {
        county_train$mask = !((county_train$fips %in% ny_counties) & (county_train$days_since_intrv > lag))
    } else {
        county_train$mask = !(county_train$fips %in% ny_counties)
    }
}

edges = read_csv("data/fips_adjacency.csv") %>% 
  # filter(isnbr) %>%
  filter(dist <= 200) %>%
  filter(src_lab %in% valid_fips, tgt_lab %in% valid_fips)

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
  spatial_scale_fixed=spatial_scale_fixed
)
saveRDS(model_data, paste0(dir, "/model_data.rds"))

model = rstan::stan_model("stan_models/1_spatiotemporal.stan")

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
  "ar_scale"
)


# pre-fitting the variational model should take ~ 1h to convergence
fit_vb = rstan::vb(
  model, 
  data=model_data,
  adapt_engaged=FALSE,
  eta = 0.35,
  iter=iter,
  tol_rel_obj=rel_tol,
  init="0",
  pars=pars,
  output_samples=samples
)
saveRDS(fit_vb, paste0(dir, "/fit.rds"))
