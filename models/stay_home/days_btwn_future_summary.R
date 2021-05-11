library(tidyverse)
library(magrittr)
library(feather)
library(rstanarm)
library(gridExtra)

up = 10
down = -10

## Read data
model <- readRDS("./model.rds")
county_future_fit <- readRDS("./county_future_fit.rds")
source("../plot_foo.R")

## Read county_future
county_pred <- read_feather("../../data/county_future_stayhome.feather")
#dim(county_pred)
#summary(county_pred$date[county_pred$index_desc == 1])

county_pred %<>%
  filter(date <= as.Date("2020-05-05"))
#dim(county_pred)

## obtain distribution values from fit sampling
county_pred %<>% 
  mutate(
    fit_mu = apply(county_future_fit, 2, mean),
    fit_med = apply(county_future_fit, 2, quantile, probs = 0.5), # use posterior median to hand skewness
    fit_lo = apply(county_future_fit, 2, quantile, probs = 0.05),
    fit_hi = apply(county_future_fit, 2, quantile, probs = 0.95))

## modify values to obtain counterfactual
county_pred1 = county_pred %>% 
  mutate(days_btwn_stayhome_thresh = days_btwn_stayhome_thresh + up) %>%
  mutate(intrv_stayhome = as.numeric(date >= stayhome + 12 + up))

county_pred3 = county_pred %>% 
  mutate(days_btwn_stayhome_thresh = days_btwn_stayhome_thresh + down) %>%
  mutate(intrv_stayhome = as.numeric(date >= stayhome + 12 + down))

county_pred1$days_since_intrv_stayhome <- county_pred1$days_since_intrv_stayhome - up
county_pred3$days_since_intrv_stayhome <- county_pred3$days_since_intrv_stayhome - down

## get posteriors
county_future_ctr1 <- model %>% 
  posterior_predict(county_pred1, draws = 500)
county_future_ctr3 <- model %>% 
  posterior_predict(county_pred3, draws = 500)

saveRDS(county_future_ctr1, "./county_future_ctr1.rds")
saveRDS(county_future_ctr3, "./county_future_ctr3.rds")

## generate nchs summaries
county_pred %<>% 
  mutate(
    ctr1_mu = apply(county_future_ctr1, 2, mean),
    ctr1_med = apply(county_future_ctr1, 2, quantile, probs = 0.5), # use posterior median to hand skewness
    ctr1_lo = apply(county_future_ctr1, 2, quantile, probs = 0.05),
    ctr1_hi = apply(county_future_ctr1, 2, quantile, probs = 0.95),
    ctr3_mu = apply(county_future_ctr3, 2, mean),
    ctr3_med = apply(county_future_ctr3, 2, quantile, probs = 0.5), # use posterior median to hand skewness
    ctr3_lo = apply(county_future_ctr3, 2, quantile, probs = 0.05),
    ctr3_hi = apply(county_future_ctr3, 2, quantile, probs = 0.95))

for(c in 1:6) {
    fips_ <- county_pred %>% 
      distinct(fips, nchs, pop) %>% 
      filter(nchs == c) %>% 
      arrange(desc(pop)) %>% 
      pull(fips)

    name_ <- county_pred %>%
      filter(fips %in% fips_) %>% 
      distinct(fips, state, county) %>%
      mutate(name = paste(state, county))

  county_plots <- lapply(fips_, 
                         function(x) county_pred %>% 
                           filter(fips == x) %>% 
                           gg_days_btwn_sampling(
                           name = name_$name[name_$fips == x], 
                           up = up, 
                           down = down, 
                           lag = 12))
  county_plots <- marrangeGrob(county_plots, 
                               nrow = 6, ncol = 2, 
                               left = "", top = "")
  ggsave(paste("./days_btwn_future_summary/", 
               "sampling_nchs_", c, ".pdf", sep = ""), 
         county_plots, width = 15, height = 25, units = "cm")
}

## generate cumulative effect summary
county_fit_effect <- matrix(nrow = 500, ncol = 0)
county_ctr1_effect <- matrix(nrow = 500, ncol = 0)
county_ctr3_effect <- matrix(nrow = 500, ncol = 0)
for(f in unique(county_pred$fips)){
  county_idx <- which(county_pred$fips == f)
  fit <- county_future_fit[, county_idx]
  fit <- t(apply(fit, 1, cumsum))
  ctr1 <- county_future_ctr1[, county_idx]
  ctr1 <- t(apply(ctr1, 1, cumsum))
  ctr3 <- county_future_ctr3[, county_idx]
  ctr3 <- t(apply(ctr3, 1, cumsum))
  
  county_fit_effect <- cbind(county_fit_effect, fit)
  county_ctr1_effect <- cbind(county_ctr1_effect, ctr1)
  county_ctr3_effect <- cbind(county_ctr3_effect, ctr3) 
}

county_pred %<>% 
  group_by(fips) %>% 
  mutate(y_eff = cumsum(y)) %>%
  ungroup()

county_pred %<>% 
  mutate(
    fit_mu_eff = apply(county_fit_effect, 2, mean),
    fit_med_eff = apply(county_fit_effect, 2, quantile, probs = 0.5), # use posterior median to hand skewness
    fit_lo_eff = apply(county_fit_effect, 2, quantile, probs = 0.05),
    fit_hi_eff = apply(county_fit_effect, 2, quantile, probs = 0.95))

county_pred %<>% 
  mutate(
    ctr1_mu_eff = apply(county_ctr1_effect, 2, mean),
    ctr1_med_eff = apply(county_ctr1_effect, 2, quantile, probs = 0.5), # use posterior median to hand skewness
    ctr1_lo_eff = apply(county_ctr1_effect, 2, quantile, probs = 0.05),
    ctr1_hi_eff = apply(county_ctr1_effect, 2, quantile, probs = 0.95),
    ctr3_mu_eff = apply(county_ctr3_effect, 2, mean),
    ctr3_med_eff = apply(county_ctr3_effect, 2, quantile, probs = 0.5), # use posterior median to hand skewness
    ctr3_lo_eff = apply(county_ctr3_effect, 2, quantile, probs = 0.05),
    ctr3_hi_eff = apply(county_ctr3_effect, 2, quantile, probs = 0.95))

for(c in 1:6) {
  fips_ <- county_pred %>% 
    distinct(fips, nchs, pop) %>% 
    filter(nchs == c) %>% 
    arrange(desc(pop)) %>% 
    pull(fips)
  
  name_ <- county_pred %>%
    filter(fips %in% fips_) %>% 
    distinct(fips, state, county) %>%
    mutate(name = paste(state, county))
  
  county_plots <- lapply(fips_, 
                         function(x) county_pred %>% 
                           filter(fips == x) %>% 
                           gg_days_btwn_effect(
                             name = name_$name[name_$fips == x], 
                             up = up, 
                             down = down, 
                             lag = 12))
  county_plots <- marrangeGrob(county_plots, 
                               nrow = 6, ncol = 2, 
                               left = "", top = "")
  ggsave(paste("./days_btwn_future_summary/", 
               "effect_nchs_", c, ".pdf", sep = ""), 
         county_plots, width = 15, height = 25, units = "cm")
}

# ## aggregate by nchs
# nchs_pred <- county_pred %>% 
#   group_by(nchs) %>% 
#   mutate(days_btwn_decrease_thresh = median(days_btwn_decrease_thresh)) %>% 
#   group_by(nchs, days_since_thresh, days_btwn_decrease_thresh) %>% 
#   summarise(y = log(median(1e-2 + county_pred$y / county_pred$pop  * 1e5)), 
#             fit_mu = NA, fit_med = NA, fit_lo = NA, fit_hi = NA, 
#             ctr1_mu = NA, ctr1_med = NA, ctr1_lo = NA, ctr1_hi = NA,
#             ctr3_mu = NA, ctr3_med = NA, ctr3_lo = NA, ctr3_hi = NA)

# county_log_fit = county_fit
# county_log_ctr = county_ctr
# for (r in 1:500) {
#   county_log_fit[r, ] = log(1e-2 + county_log_fit[r, ] / county_pred$pop * 1e5)
#   county_log_ctr[r, ] = log(1e-2 + county_log_ctr[r, ] / county_pred$pop * 1e5)
# }

# for(n in unique(nchs_pred$nchs)){
#   days <- nchs_pred$days_since_thresh[which(nchs_pred$nchs == n)]
#   for(d in days) {
#     county_idx <- which(county_pred$nchs == n & 
#                           county_pred$days_since_thresh == d)
#     if(length(county_idx) > 1) {
#       fit <- apply(county_log_fit[, county_idx], 1 , mean)
#       ctr <- apply(county_log_ctr[, county_idx], 1 , mean)
#     } else {
#       fit <- county_log_fit[, county_idx]
#       ctr <- county_log_ctr1[, county_idx]
#     }
#     nchs_idx <- which(nchs_pred$nchs == n & 
#                         nchs_pred$days_since_thresh == d)
    
#     nchs_pred$fit_mu[nchs_idx] <- mean(fit)
#     nchs_pred$fit_med[nchs_idx] <- quantile(fit, 0.5)
#     nchs_pred$fit_lo[nchs_idx] <- quantile(fit, 0.05)
#     nchs_pred$fit_hi[nchs_idx] <- quantile(fit, 0.95)
    
#     nchs_pred$ctr_mu[nchs_idx] <- mean(ctr)
#     nchs_pred$ctr_med[nchs_idx] <- quantile(ctr, 0.5)
#     nchs_pred$ctr_lo[nchs_idx] <- quantile(ctr, 0.05)
#     nchs_pred$ctr_hi[nchs_idx] <- quantile(ctr, 0.95)
#   }
# }

# county_plots <- lapply(1:6, 
#                        function(x) nchs_pred %>% 
#                          filter(nchs == x)%>% 
#                          gg_intrv_agg_sampling(name = x, 
#                                            intrv_name = "decrease", 
#                                            lag = 5))
# county_plots <- marrangeGrob(county_plots, 
#                              nrow = 6, ncol = 2, 
#                              left = "", top = "")
# ggsave(paste("./intervention_0_summary/", 
#              "nchs_sampling.pdf", sep = ""), 
#        county_plots, width = 15, height = 25, units = "cm")
