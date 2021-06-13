library(feather)
library(tidyverse)
library(lubridate)
library(ape)
library(argparse)
library(gridExtra)
library(rstan)
library(cowplot)

source("./utils.R")


parser = ArgumentParser()
parser$add_argument("--dir", type="character", default="100k_ar1", 
    help="Directory where results where saved")


args = parser$parse_args()
for (i in 1:length(args)) {
    varname = names(args)[i]
    value = args[[i]]
    if (varname == "dir")
        value = paste0("results/", value)
    assign(varname, value, envir = .GlobalEnv)
}
dir.create(paste0(dir, "/residuals"), recursive=TRUE)

fit = read_rds(sprintf("%s/fit.rds", dir))
model_data = read_rds(sprintf("%s/model_data.rds", dir))

# compute and standardize residuals
y = model_data$df$y
pred = exp(rstan::extract(fit, pars="log_rate")[[1]])

# form 1 of standardization
# pred = exp(predicted$log_yhat)
yhat_ = apply(pred, 2, mean)
residuals = - (yhat_ - y)

y_ = matrix(rep(y, each=nrow(pred)), nrow=nrow(pred))
deviance = 2 * y_ * log((y_ + 1e-1) / (pred + 1e-1)) + 2 * (1 + y_) * log((1 + pred) / (1 + y_))
residuals_rel = apply(sign(y_ - pred) * sqrt(abs(deviance)), 2, mean)

# 1. Autocorrelation histogram ----------

max_cum_deaths = model_data$df %>% 
  group_by(fips) %>% 
  summarize(mdp = max(cum_deaths), .groups="drop") %>% 
  arrange(desc(mdp))
fips = max_cum_deaths$fips[1:50]
acfs = numeric(0)
for (i in 1:length(fips)) {
  ix = which(model_data$df$fips == fips[i])
  rrel = residuals_rel[ix] + 0.01 * rnorm(length(ix))
  acfs[i] = acf(rrel, plot = FALSE)$acf[2]
}

{
  png(sprintf("%s/residuals/acf.png", dir))
  hist(acfs,  xlim=c(-1, 1), breaks=10)
  abline(v=mean(acfs), lty=2, col="red")
  print(paste("Mean autocorrelation in largest 50 counties: ", mean(acfs)))
  dev.off()
}

#Load model residuals
res_df = tibble(
  fips=model_data$df$fips,
  res=residuals,
  y=model_data$df$y,
  yhat=yhat_,
  res_rel=residuals_rel,  # (mu - y) / mu
  t=model_data$df$days_since_thresh
)

res_plots <- list()
for(t_ in 7:22) {
  res_df_t = res_df %>% 
  filter(t == t_)
  
  mn = mean(res_df_t$res)
  lines = tibble(
    v=c(quantile(res_df_t$res, .05), 
        quantile(res_df_t$res, .95)),
    type=c("5%", "95%"))
  
  res_plots[[as.character(t_)]] <- ggplot(res_df_t) + 
    geom_histogram(aes(x=res), bins=60) + 
    geom_vline(aes(xintercept=mn), color="blue") +
    geom_vline(data=lines, aes(xintercept=v, lty=type), color="red") +
    theme_bw() +
    labs(title=sprintf("Residuals for epidemic time %d", t_),
         color="")
}

res_plots <- marrangeGrob(res_plots, 
                          nrow = 8, ncol = 2, 
                          left = "", top = "")

ggsave(sprintf('%s/residuals/res_plots.pdf', dir), res_plots, width = 8 height = 8 units = "in")


res_rel_plots <- list()
for(t_ in 7:22) {
  res_df_t = res_df %>% 
  filter(t == t_) %>%
  filter(abs(res_rel) > 1e-2)
  
  mn = mean(res_df_t$res_rel)
  lines = tibble(
    v=c(quantile(res_df_t$res_rel, .05), 
        quantile(res_df_t$res_rel, .95)),
    type=c("5%", "95%"))
  
  res_rel_plots[[as.character(t_)]] <- ggplot(res_df_t) + 
    geom_histogram(aes(x=res_rel), bins=20) + 
    geom_vline(aes(xintercept=mn), color="blue") +
    geom_vline(data=lines, aes(xintercept=v, lty=type), color="red") +
    theme_bw() +
    labs(title=sprintf("Residuals for epidemic time %d", t_),
         color="", linetype="", x="residuals", y="") +
    theme(title=element_text(size=8))
}

res_rel_plots <- marrangeGrob(res_rel_plots, 
                          nrow = 8, ncol = 2, 
                          left = "", top = "")

ggsave(sprintf('%s/residuals/res_rel_plots.pdf', dir),res_rel_plots, width = 8 height = 8 units = "in")


distmat_raw = read_rds("data/OD.rds")
distmat_raw = distmat_raw + t(distmat_raw)

adjacency = read_csv("data/edges_with_knn2_filt100.csv")


# moran's I

moransI_list <- list()  
for(t_ in 7:22)
{
  res_df_t = res_df %>% 
    filter(t == t_)
  x = res_df_t$res
    #q05 = quantile(res_df_t$res, .05)
  insample = x > -Inf
  x = x[insample]
  fips = res_df_t$fips[insample]
  wts = 1 / distmat_raw[fips, fips]
  diag(wts) = 0

  # From package
  moransI_list[[t_]] <- as.data.frame(Moran.I(x, wts, scaled=TRUE))  # from package
  moransI_list[[t_]]$t = t_
}
moransI = bind_rows(moransI_list)

moransI %>% 
  pivot_longer(-t) %>% 
  ggplot() + 
  geom_line(aes(x=t, y=value)) +
  facet_wrap(~name, scales = "free", ncol=1) +
  theme_minimal_hgrid()
ggsave(sprintf('%s/residuals/moransI_distance.pdf', dir), width = 8 height = 8 units = "in")


# ---------------
#summary(res_df$t)
moransI_list <- list()  
for(t_ in 7:22)
{
  res_df_t = res_df %>% 
    filter(t == t_)
  x = res_df_t$res_rel

  insample = x > -Inf #q05 = quantile(res_df_t$res, .05)
  x = x[insample]
  fips = res_df_t$fips[insample]
  wts = 1 / distmat_raw[fips, fips]
  diag(wts) = 0

  # From package
  moransI_list[[t_]] <- as.data.frame(Moran.I(x, wts, scaled=TRUE))  # from package
  moransI_list[[t_]]$t = t_
}
moransI = bind_rows(moransI_list)

moransI %>% 
  pivot_longer(-t) %>% 
  ggplot() + 
  geom_line(aes(x=t, y=value)) +
  facet_wrap(~name, scales = "free", ncol=1) +
  theme_minimal_hgrid()
ggsave(sprintf('%s/residuals/moransI_distance_relative.pdf', dir), width = 8 height = 8 units = "in")


# Adjacency
#summary(res_df$t)
moransI_list <- list()  
for(t_ in 7:22)
{
  res_df_t = res_df %>% 
    filter(t == t_)
  x = res_df_t$res
    
  insample = x > -Inf #q05 = quantile(res_df_t$res, .05)
  x = x[insample]
  fips = res_df_t$fips[insample]
  
  wts = matrix(0, length(fips), length(fips), dimnames = list(fips, fips))
  adj_t = adjacency %>% 
    filter(src_lab %in% fips, tgt_lab %in% fips)
  num_edges = nrow(adj_t)
  for(i in 1:num_edges) {
    s = adj_t$src_lab[i]
    t = adj_t$tgt_lab[i]
    wts[s, t] = wts[t, s] = 1
  }

  # From package
  moransI_list[[t_]] <- as.data.frame(Moran.I(x, wts, scaled=TRUE))  # from package
  moransI_list[[t_]]$t = t_
}
moransI = bind_rows(moransI_list)

moransI %>% 
  pivot_longer(-t) %>% 
  ggplot() + 
  geom_line(aes(x=t, y=value)) +
  facet_wrap(~name, scales = "free", ncol=1) +
  theme_minimal_hgrid()
ggsave(sprintf('%s/residuals/moransI_adjacency.pdf', dir), width = 8 height = 8 units = "in")


moransI_list <- list()  
for(t_ in 7:22)
{
  res_df_t = res_df %>% 
    filter(t == t_)
  x = res_df_t$res_rel
    
  insample = x > -Inf #q05 = quantile(res_df_t$res, .05)
  x = x[insample]
  fips = res_df_t$fips[insample]
  
  wts = matrix(0, length(fips), length(fips), dimnames = list(fips, fips))
  adj_t = adjacency %>% 
    filter(src_lab %in% fips, tgt_lab %in% fips)
  num_edges = nrow(adj_t)
  for(i in 1:num_edges) {
    s = adj_t$src_lab[i]
    t = adj_t$tgt_lab[i]
    wts[s, t] = wts[t, s] = 1
  }

  # From package
  moransI_list[[t_]] <- as.data.frame(Moran.I(x, wts, scaled=TRUE))  # from package
  moransI_list[[t_]]$t = t_
}
moransI = bind_rows(moransI_list)

moransI %>% 
  pivot_longer(-t) %>% 
  ggplot() + 
  geom_line(aes(x=t, y=value)) +
  facet_wrap(~name, scales = "free", ncol=1) +
  theme_minimal_hgrid()
ggsave(sprintf('%s/residuals/moransI_adjacency_relative.pdf', dir), width = 8 height = 8 units = "in")
