library(tidyverse)
library(cowplot)
library(feather)
library(igraph)
library(lubridate)
library(yaml)
options(mc.cores = 4)
source("./utils.R")

what = "stayhome"

fit_full = read_rds(sprintf("./results/%s/vb/full_model_no_temporal/fit.rds", what))
fit_no_ny = read_rds(sprintf("./results/%s/vb/no_temporal_no_ny/fit.rds", what))

model_data_full = read_rds(sprintf("./results/%s/vb/full_model/model_data.rds", what))
model_data_no_ny = read_rds(sprintf("./results/%s/vb/no_temporal_no_ny/model_data.rds", what))

parnames = c(
  "nchs_pre", "nchs_post", "beta_covars_pre",
  "beta_covars_post", "beta_covars_pre_inter",
  "state_eff"
)

pars_full = rstan::extract(fit_full)
pars_no_ny = rstan::extract(fit_no_ny)

iqrfun = function(x) quantile(x, 0.75) - quantile(x, 0.25)

dfs = list()
k = 0

new_names = c(
  beta_covars_post="post-interaction",
  beta_covars_pre="pre-interaction",
  beta_covars_pre_inter="pre-interaction"
)
# p = "beta_covars_post"
for (p in c("beta_covars_post", "beta_covars_pre", "beta_covars_pre_inter")) {
  p_full = apply(pars_full[[p]], c(2, 3), median)
  p_no_ny = apply(pars_no_ny[[p]], c(2, 3), median)
  p_full_iqr = apply(pars_full[[p]], c(2,3), iqrfun)
  p_no_ny_iqr = apply(pars_no_ny[[p]], c(2, 3), iqrfun)
  
  k = k + 1
  M = length(as.numeric(p_full))
  timing_ind = numeric(M)
  R = nrow(p_full)
  C = ncol(p_full)
  if (k == 1)
    for (c in 1:C)
      timing_ind[1 + R * (c - 1)] = 1
  dfs[[k]] = tibble(
    parname = ifelse(timing_ind == 0, new_names[p], "intervention timing"),
    full_median = as.numeric(p_full),
    full_iqr = as.numeric(p_full_iqr),
    no_ny_median = as.numeric(p_no_ny),
    no_ny_iqr = as.numeric(p_no_ny_iqr),
    timing_ind = timing_ind
  )
}
dfs = bind_rows(dfs)

ggplot(dfs) +
  geom_point(
    aes(
      x=full_median,
      y=no_ny_median,
      shape=parname,
      size=timing_ind,
      color=parname
    ),
    alpha=0.5
  ) +
  geom_abline(slope=1, intercept=0, col="black", lty=2) + 
  scale_color_manual(values=c("red", "black", "black")) +
  scale_size_continuous(range = c(2,4)) +
  guides(size=FALSE) +
  labs(shape="Parameter", color="Parameter", x="NYC included", y="NYC removed") +
  theme_bw()
ggsave("review/new_york/params_comparison.pdf", width=5, height=3, units="in")


dat_full = read_csv(sprintf("results/%s/vb/full_model_no_temporal/counterfactuals/plotdata_early_late_action_totals_no_ny.csv", what))
dat_no_ny = read_csv(sprintf("results/%s/vb/no_temporal_no_ny/counterfactuals/plotdata_early_late_action_totals_no_ny.csv", what))

plotdata = bind_rows(
  mutate(dat_full, src="NYC included"),
  mutate(dat_no_ny, src="NYC removed")
) %>% 
  filter(date > lubridate::ymd("20200325"))

uplim = 30000
dolim = -30000

ggplot(plotdata) +
  geom_line(aes(x=date, y=pmax(dolim, pmin(med, uplim)), linetype=src, color=action)) +
  geom_ribbon(aes(x=date, ymax=pmax(dolim, pmin(hi, uplim)), linetype=src, ymin=pmax(dolim, pmin(lo, uplim)), fill=action), alpha=0.15) +
  theme_minimal_hgrid() +
  scale_color_manual(values=c("#009E73", "#D55E00")) +
  scale_fill_manual(values=c("#009E73", "#D55E00")) +
  theme(axis.text.y = element_text(size=8)) +
  labs(fill="Timing", color="Timing", y = "Excess/Averted Deaths", color="", linetype="Dataset") + 
  theme(axis.title.x = element_blank()) +
  guides()

ggsave("review/new_york/counterfactual_comparison.pdf", width=6, height=3, units="in")
