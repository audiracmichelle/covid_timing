library(tidyverse)
library(magrittr)
library(feather)
library(gridExtra)
library(cowplot)
library(reshape2)
library(rstan)
library(argparse)
options(mc.cores = 4)
source("./utils.R")


parser = ArgumentParser()
parser$add_argument("--dir", type="character", default="100k_ar1", 
    help="Directory where results where saved")
parser$add_argument("--shift_days", type="integer", default=10, 
    help="How much to shift up and down.")
parser$add_argument("--exclude_ny", action="store_true", default=FALSE, 
    help="Remove NYC counties for counterfactual plots.")


args = parser$parse_args()
for (i in 1:length(args)) {
    varname = names(args)[i]
    value = args[[i]]
    if (varname == "dir")
        value = paste0("results/", value)
    assign(varname, value, envir = .GlobalEnv)
}
dir.create(paste0(dir, "/counterfactuals"), recursive=TRUE)
# dir = "../results/decrease/vb/full_model"
suffix = ifelse(exclude_ny, "no_ny", "")
fit = read_rds(sprintf("%s/fit.rds", dir))
model_data = read_rds(sprintf("%s/model_data.rds", dir))

ny_counties = c(
  "36081",  # queens
  "36085",  # richmond
  "36061",  # new york
  "36047",  # kings
  "36005"  # bronx
)

bent_cable=model_data$bent_cable
rand_eff=TRUE
states=TRUE
spatial=model_data$spatial
temporal = model_data$temporal

up = shift_days
down = -shift_days
intervention = model_data$type

evaldata = read_feather("data/county_train_.feather") %>%   # 1454 counties
  filter(fips %in% levels(model_data$df$fips_f)) %>%   # 1021 counties
  filter(days_since_thresh <= 30) %>% 
  mutate(interv_type = .env$intervention) %>%
  mutate(days_since_intrv = if_else(interv_type == "decrease", days_since_intrv_decrease, days_since_intrv_stayhome)) %>%
  mutate(days_btwn = if_else(interv_type == "decrease", days_btwn_decrease_thresh, days_btwn_stayhome_thresh)) %>%
  group_by(fips) %>% 
  filter(
    max(days_since_thresh) >= 7,
    max(y) >= 1  # for evaluation max(cum_deaths >= 1) is not a good idea cause there are tons of zeros
  ) %>%  
  ungroup()

if (exclude_ny)
  evaldata = filter(evaldata, !(fips %in% ny_counties))

# evaldata = model_data$df

county_fit_posterior = posterior_predict(
  fit,
  model_data,
  new_df=evaldata,
  rand_eff=rand_eff,
  bent_cable=bent_cable,
  states=states,
  spatial=spatial,
  temporal=FALSE
)
yhat_no_rand_eff = exp(county_fit_posterior$log_yhat) # / exp( county_fit_posterior$rand_eff)
county_fit = yhat_no_rand_eff %>% 
  melt() %>% 
  `names<-`(c("sim", "row", "yhat")) %>% 
  left_join(
    select(
      mutate(evaldata, row=1:n()),
      row, nchs, days_since_thresh, days_since_intrv, pop
    )
  ) %>% 
  mutate(yhat = 1e6 * yhat/pop)

county_fit_posterior_down = posterior_predict(
  fit,
  model_data,
  new_df=evaldata,
  rand_eff=rand_eff,
  states=states,
  spatial=spatial,
  bent_cable=bent_cable,
  shift_timing=-shift_days,
  temporal=FALSE
)
yhat_no_rand_eff = exp(county_fit_posterior_down$log_yhat) # / exp(county_fit_posterior$rand_eff)
county_fit_down = yhat_no_rand_eff %>% 
  melt() %>% 
  `names<-`(c("sim", "row", "yhat")) %>% 
  left_join(
    select(
      mutate(evaldata, row=1:n(), days_since_intrv=days_since_intrv - down),
      row, nchs, days_since_thresh, days_since_intrv, pop
    )
  )  %>% 
  mutate(yhat =1e6 *  yhat/pop)

county_fit_posterior_up = posterior_predict(
  fit,
  model_data,
  new_df=evaldata,
  rand_eff=rand_eff,
  states=states,
  spatial=spatial,
  bent_cable=bent_cable,
  temporal=FALSE,
  shift_timing=shift_days
)
yhat_no_rand_eff = exp(county_fit_posterior_up$log_yhat)
county_fit_up = yhat_no_rand_eff %>% 
  melt() %>% 
  `names<-`(c("sim", "row", "yhat")) %>% 
  left_join(
    select(
      mutate(evaldata, row=1:n(), days_since_intrv=days_since_intrv - up),
      row, nchs, days_since_thresh, days_since_intrv, pop
    )
  )  %>% 
  mutate(yhat = 1e6 * yhat/pop)

county = county_fit %>%
  group_by(nchs, sim, days_since_thresh) %>%
  summarize(yhat = median(yhat), .groups="drop") %>% 
  mutate(log_yhat=log(yhat)) %>% 
  group_by(nchs, days_since_thresh) %>% 
  summarize(
    fit_mu = mean(yhat),
    fit_med = median(yhat), # use posterior median to hand skewness
    fit_lo = quantile(yhat, 0.05),
    fit_hi = quantile(yhat, 0.95),
    .groups="drop"
  ) %>% 
  mutate(type="actual")

county_up = county_fit_up %>%
  group_by(nchs, sim, days_since_thresh) %>%
  summarize(yhat = median(yhat), .groups="drop") %>% 
  mutate(log_yhat=log(yhat)) %>% 
  group_by(nchs, days_since_thresh) %>% 
  summarize(
    fit_mu = mean(yhat),
    fit_med = median(yhat), # use posterior median to hand skewness
    fit_lo = quantile(yhat, 0.05),
    fit_hi = quantile(yhat, 0.95),
    .groups="drop"
  ) %>% 
  mutate(type="late")

county_down = county_fit_down  %>%
  group_by(nchs, sim, days_since_thresh) %>%
  summarize(yhat = median(yhat), .groups="drop") %>% 
  mutate(log_yhat=log(yhat)) %>% 
  group_by(nchs, days_since_thresh) %>% 
  summarize(
    fit_mu = mean(yhat),
    fit_med = median(yhat), # use posterior median to hand skewness
    fit_lo = quantile(yhat, 0.05),
    fit_hi = quantile(yhat, 0.95),
    .groups="drop"
  ) %>%
  mutate(type="early")

  plotdata = bind_rows(
  county,
  county_up,
  county_down
) %>% 
  filter(days_since_thresh <= 30)

uplim = 30
dolim = 0.01
lagat = 14

ggplot(plotdata) +
  geom_line(aes(x=days_since_thresh, y=pmax(dolim, pmin(fit_med, uplim)), color=type), lwd=1) +
#   geom_vline(aes(xintercept=days_btwn_thresh + lagat), data=interv_day, lty=2) +
  geom_ribbon(aes(x=days_since_thresh, ymax=pmax(dolim, pmin(fit_hi, uplim)), ymin=pmax(dolim, pmin(fit_lo, uplim)), fill=type), alpha=0.4) +
  scale_y_log10() +
  # facet_wrap(~ paste("NCHS", nchs)) +
  facet_wrap(~ paste("NCHS", nchs), scales = "free_y") +
    guides(color=FALSE) +
  # scale_y_continuous(
  #   limits=c(-5.5, 2.2),
  #   breaks=c(-1, 0, 1),
  #   labels=c("0.1", "1", "10")
  # ) +
  theme_minimal_hgrid() +
  scale_color_manual(values=c("#0072B2", "#009E73", "#D55E00")) +
  scale_fill_manual(values=c("#0072B2", "#009E73", "#D55E00")) +
  theme(legend.position = "top") +
  labs(fill="", y = "Deaths per 1 million", x="Days since threshold deaths")

ggsave(sprintf("%s/counterfactuals/counterfactuals_%s.pdf", dir, suffix), width=8,height=4,units="in")


counterf_eval_data = model_data$df
# if (exclude_ny)
#   counterf_eval_data = filter(counterf_eval_data, !(fips %in% ny_counties))

# temporal = TRUE  # override for now
# diffs in diffs counterfactuals
y_samples_actual = posterior_predict(
  fit,
  model_data,
  new_df=counterf_eval_data,
  rand_eff=rand_eff,
  states=states,
  spatial=spatial,
  bent_cable=bent_cable,
  temporal=temporal
)$y_samples %>%
  melt() %>% 
  `names<-`(c("sim", "row", "actual")) %>% 
  left_join(mutate(ungroup(counterf_eval_data), row=1:n())) %>%
  # group_by(sim, fips) %>%
  # arrange(sim, fips, date) %>%
  # mutate(actual=cumsum(actual)) %>%
  # ungroup() %>%
  na.omit() %>% 
  filter(!(fips %in% ny_counties))


y_samples_down = posterior_predict(
  fit,
  model_data,
  new_df=counterf_eval_data,
  rand_eff=rand_eff,
  states=states,
  spatial=spatial,
  bent_cable=bent_cable,
  temporal=temporal,
  shift_timing=-shift_days
)$y_samples %>%
  melt() %>% 
  `names<-`(c("sim", "row", "down")) %>% 
  left_join(mutate(ungroup(counterf_eval_data), row=1:n())) %>%
  # group_by(sim, fips) %>%
  # arrange(sim, fips, date) %>%
  # mutate(down=cumsum(down)) %>%
  # ungroup()  %>%
  na.omit() %>% 
  filter(!(fips %in% ny_counties))



y_samples_up = posterior_predict(
  fit,
  model_data,
  new_df=counterf_eval_data,
  rand_eff=rand_eff,
  states=states,
  spatial=spatial,
  bent_cable=bent_cable,
  temporal=temporal,
  shift_timing=shift_days
)$y_samples %>%
  melt() %>% 
  `names<-`(c("sim", "row", "up")) %>% 
  left_join(mutate(ungroup(counterf_eval_data), row=1:n())) %>%
  # group_by(sim, fips) %>%
  # arrange(sim, fips, date) %>%
  # mutate(up=cumsum(up)) %>%
  # ungroup() %>%
  na.omit() %>% 
  filter(!(fips %in% ny_counties))



county_late_action = y_samples_actual %>%
  left_join(y_samples_down) %>%
  left_join(y_samples_up) %>% 
  mutate(late_action=up - actual) %>%
  group_by(nchs, sim, date) %>%
  summarize(late_action=sum(late_action), .groups="drop") %>%
  group_by(nchs, sim) %>%
  mutate(late_action=cumsum(late_action)) %>%
  ungroup() %>%
  group_by(nchs, date) %>%
  summarize(
    mu = mean(late_action),
    med = median(late_action), # use posterior median to hand skewness
    lo = quantile(late_action, 0.05),
    hi = quantile(late_action, 0.95),
    .groups="drop"
  ) %>%
  mutate(action="10 days later")

county_early_action = y_samples_actual %>%
  left_join(y_samples_down) %>%
  left_join(y_samples_up) %>% 
  mutate(early_action=down - actual) %>%
  group_by(nchs, sim, date) %>%
  summarize(early_action=sum(early_action), .groups="drop") %>%
  group_by(nchs, sim) %>%
  mutate(early_action=cumsum(early_action)) %>%
  ungroup() %>%
  group_by(nchs, date) %>%
  summarize(
    mu = mean(early_action),
    med = median(early_action), # use posterior median to hand skewness
    lo = quantile(early_action, 0.05),
    hi = quantile(early_action, 0.95),
    .groups="drop"
  ) %>%
  mutate(action="10 days earlier")

 plotdata = bind_rows(county_late_action, county_early_action) %>%
  filter(date >= lubridate::ymd("20200315")) %>%
  mutate(nchs=paste("NCHS", nchs))

  dolim=-1e6
  uplim =1e6

weird = scales::trans_new(
  "signed_log",
  transform=function(x) sign(x)*log(1 + abs(x)),
  inverse=function(x) sign(x)*(exp(abs(x)) - 1)
)

# weird = scales::trans_new(
#   "signed_sqrt",
#   transform=function(x) sign(x)*sqrt(abs(x)),
#   inverse=function(x) sign(x)*(abs(x))^2
# )

interv_day = model_data$df %>%
  # filter(fips %in% valid_fips) %>%   # 1021 counties
  # filter(days_since_thresh <= 30) %>%
  select(nchs, days_btwn_decrease_thresh) %>%
  group_by(nchs) %>%
  summarize(
    days_btwn_decrease_thresh=mean(days_btwn_decrease_thresh, na.rm=TRUE),
    .groups="drop")  %>%
  mutate(nchs=paste("NCHS", nchs))
  

print(max(plotdata$med))

lagat = 14
ggplot(plotdata) +
  geom_line(aes(x=date, y=pmax(dolim, pmin(med, uplim)), color=action)) +
  geom_ribbon(aes(x=date, ymax=pmax(dolim, pmin(hi, uplim)), ymin=pmax(dolim, pmin(lo, uplim)), fill=action), alpha=0.4) +
  facet_wrap(~nchs, ncol=3, scale="free_y") +
  theme_minimal_hgrid() +
  # geom_vline(aes(xintercept=days_btwn_decrease_thresh + lagat), data=interv_day, lty=2) +
  # scale_y_continuous(trans=weird, n.breaks=8) + # , breaks=c(-15000, -10000, -5000, -1000, -500, -100, 0, 100, 500 1000, 5000, 10000, 15000)) +
  scale_color_manual(values=c("#009E73", "#D55E00")) +
  scale_fill_manual(values=c("#009E73", "#D55E00")) +
  theme(legend.position = "top", axis.text.y = element_text(size=8)) +
  labs(fill="", y = "Excess/Averted Deaths", color="") + 
  theme(axis.title.x = element_blank(), legend.position="top") +
  guides(fill=FALSE)


ggsave(sprintf("%s/counterfactuals/early_late_action_%s.pdf", dir, suffix), width=8,height=4,units="in")
write_csv(plotdata, sprintf("%s/counterfactuals/plotdata_early_late_action_%s.csv", dir, suffix))


# same but aggregated

county_early_action_tot = y_samples_actual %>%
  left_join(y_samples_down) %>%
  left_join(y_samples_up) %>% 
  mutate(early_action=down - actual) %>%
  group_by(sim, date) %>%
  summarize(early_action=sum(early_action), .groups="drop") %>%
  group_by(sim) %>%
  mutate(early_action=cumsum(early_action)) %>%
  ungroup() %>%
  group_by(date) %>%
  summarize(
    mu = mean(early_action),
    med = median(early_action), # use posterior median to hand skewness
    lo = quantile(early_action, 0.05),
    hi = quantile(early_action, 0.95),
    .groups="drop"
  ) %>%
  mutate(action="10 days earlier")

  county_late_action_tot = y_samples_actual %>%
    left_join(y_samples_down) %>%
    left_join(y_samples_up) %>% 
    mutate(late_action=up - actual) %>%
    group_by(sim, date) %>%
    summarize(late_action=sum(late_action), .groups="drop") %>%
    group_by(sim) %>%
    mutate(late_action=cumsum(late_action)) %>%
    ungroup() %>%
    group_by(date) %>%
    summarize(
      mu = mean(late_action),
      med = median(late_action), # use posterior median to hand skewness
      lo = quantile(late_action, 0.05),
      hi = quantile(late_action, 0.95),
      .groups="drop"
    ) %>%
    mutate(action="10 days later")

    
 plotdata = bind_rows(county_late_action_tot, county_early_action_tot) %>%
  filter(date >= lubridate::ymd("20200315"))

ggplot(plotdata) +
  geom_line(aes(x=date, y=pmax(dolim, pmin(med, uplim)), color=action)) +
  geom_ribbon(aes(x=date, ymax=pmax(dolim, pmin(hi, uplim)), ymin=pmax(dolim, pmin(lo, uplim)), fill=action), alpha=0.4) +
  # facet_wrap(~nchs, ncol=3, scale="free_y") +
  theme_minimal_hgrid() +
  # geom_vline(aes(xintercept=days_btwn_decrease_thresh + lagat), data=interv_day, lty=2) +
  # scale_y_continuous(trans=weird, n.breaks=8) + # , breaks=c(-15000, -10000, -5000, -1000, -500, -100, 0, 100, 500 1000, 5000, 10000, 15000)) +
  scale_color_manual(values=c("#009E73", "#D55E00")) +
  scale_fill_manual(values=c("#009E73", "#D55E00")) +
  theme(legend.position = "top", axis.text.y = element_text(size=8)) +
  labs(fill="", y = "Excess/Averted Deaths", color="") + 
  theme(axis.title.x = element_blank(), legend.position="top") +
  guides(fill=FALSE)

ggsave(sprintf("%s/counterfactuals/early_late_action_totals.pdf", dir), width=8,height=4,units="in")
write_csv(plotdata, sprintf("%s/counterfactuals/plotdata_early_late_action_totals_%s.csv", dir, suffix))
