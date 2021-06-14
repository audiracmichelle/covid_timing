library(tidyverse)
library(rstan)
library(cowplot)
library(argparse)
options(mc.cores = 4)
source("./utils.R")

parser = ArgumentParser()
parser$add_argument("--dir", type="character", default="100k_ar1", 
    help="Directory where results where saved")
# parser$add_argument("--bent_cable", action="store_true", default=FALSE,
#     help="Forces the fit to the model with random lag")

args = parser$parse_args()
for (i in 1:length(args)) {
    varname = names(args)[i]
    value = args[[i]]
    if (varname == "dir")
        value = paste0("results/", value)
    assign(varname, value, envir = .GlobalEnv)
}
dir.create(paste0(dir, "/summaries"), recursive=TRUE)
# dir = "results/decrease/vb/no_temporal_no_cities_post"

fit = read_rds(sprintf("%s/fit.rds", dir))
model_data = read_rds(sprintf("%s/model_data.rds", dir))

ysamples = posterior_predict(
    fit,
    model_data,
    bent_cable=model_data$bent_cable,
    temporal=model_data$temporal,
    spatial=model_data$spatial
)$y_samples
postpred = posterior_predict(
    fit,
    model_data,
    bent_cable=model_data$bent_cable,
    temporal=FALSE,
    spatial=model_data$spatial
)
county_lp_var = exp(postpred$log_yhat)
prev_trend = exp(postpred$pre_term)
# county_lp_var = exp(rstan::extract(fit, pars="log_rate")[[1]])
# county_lp_var = exp(posterior_predict(fit, model_data, spatial=TRUE, bent_cable=bent_cable, temporal=FALSE)$log_yhat)
# county_lp_var = ysamples

county_train = model_data$df
fiplist = c(
    "Los Angeles"="06037",
    "Kings (WA)"="53033",
    "Wayne (Michigan)"="26163",
    "New York"="36061",
    "Queens"="36081"
)
# f1 = "06037"
use_mask = model_data$use_mask

for (i in 1:length(fiplist)) {
    f1 = fiplist[i]
    ix = which(model_data$df$fips == f1)
    if (length(ix) == 0)
        next
    pop = 1
    yi = model_data$y[ix] / pop
    predi = ysamples[, ix] / pop
    predi2 = county_lp_var[ ,ix] / pop
    yhati = apply(predi, 2, mean)
    yhati2 = apply(predi2, 2, mean)
    previ = apply(prev_trend[,ix]/ pop, 2, mean)
    yhati_95 = apply(predi, 2, quantile, .95)
    yhati_05 = apply(predi, 2, quantile, .05)
    yhati_95_trend = apply(predi2, 2, quantile, .95)
    yhati_05_trend = apply(predi2, 2, quantile, .05)
    {
        png(sprintf("%s/summaries/curve_%s.png", dir, f1))
        plot(yi, ylim=c(min(yhati_05), max(yhati_95)))
        lines(yhati, col="red", lty=2)
        lines(yhati2, col="red", lty=1)
        lines(yhati_95, col="blue", lty=2)
        lines(yhati_05, col="blue", lty=2)
        # dbtwn = model_data$df[ix, 1]
        # abline(v=dbtwn + 15, lty=3, col="gray")
        title(sprintf("FIPS %s", names(fiplist)[i]))
        dev.off()
    }

    plotdata = tibble(
        med = yhati,
        trend = yhati2, 
        y = yi,
        previ = previ,
        hi = yhati_95,
        lo = yhati_05,
        hi_trend = yhati_95_trend,
        lo_trend = yhati_05_trend,
        mask = model_data$mask[ix],
        date = model_data$df$date[ix]
    )
    ggplot(plotdata) +
        # geom_line(aes(x=date, y=Median, linetype="Y")) +
        geom_line(aes(x=date, y=trend, linetype="Pre + Post")) +
        geom_line(aes(x=date, y=previ, linetype="Pre")) +
        geom_ribbon(aes(x=date, ymax=hi, ymin=lo, fill="90% CI"), alpha=0.1) +
        # geom_ribbon(aes(x=date, ymax=hi_trend, ymin=lo_trend, fill="Trend"), alpha=0.1) +
        geom_point(aes(x=date, y=y, shape=factor(mask, levels=c(0, 1))), size=2) +
        theme_minimal_hgrid() +
        theme(axis.title.x = element_blank()) +
        labs(y = "Daily Deaths", linetype="Trend", shape=ifelse(use_mask, "Data", ""), fill="") +
        guides(fill=FALSE) +
        scale_linetype_manual(values=c(2, 1)) +
        scale_shape_manual(values=c(21, 19), labels=c("Observed", "Heldout")) #+
        # scale_fill_manual(values=c("red", "blue"), labels=c("Negative Binomial Fit (Y)", expression("Trend (N  lambda")))
    ggsave(sprintf("%s/summaries/curve_%s_pre_post.png", dir, f1), width=6, height=4, units="in")
        
}

{
    png(sprintf("%s/summaries/spatial_scale.png", dir))
    par(mfrow=c(1,3))
    hist(rstan::extract(fit, "spatial_scale[1]")[[1]])
    hist(rstan::extract(fit, "spatial_scale[2]")[[1]])
    hist(rstan::extract(fit, "spatial_scale[3]")[[1]])
    dev.off()
}


{
    png(sprintf("%s/summaries/rand_eff_scale.png", dir))
    par(mfrow=c(1,3))
    hist(rstan::extract(fit, "scale_rand_eff[1]")[[1]])
    hist(rstan::extract(fit, "scale_rand_eff[2]")[[1]])
    hist(rstan::extract(fit, "scale_rand_eff[3]")[[1]])
    dev.off()
}

{
    png(sprintf("%s/summaries/rand_eff.png", dir))
    par(mfrow=c(1,3))
    rand_eff = rstan::extract(fit, "rand_eff")[[1]]
    rand_eff = apply(rand_eff,c(2,3), median)
    hist(rand_eff[ ,1])
    hist(rand_eff[ ,2])
    hist(rand_eff[ ,3]) 
    dev.off()
}

{
    png(sprintf("%s/summaries/spatial_eff.png", dir))
    par(mfrow=c(1,3))
    spatial_eff = rstan::extract(fit, "spatial_eff")[[1]]
    spatial_eff = apply(spatial_eff,c(2,3), median)
    hist(spatial_eff[ ,1])
    hist(spatial_eff[ ,2])
    hist(spatial_eff[ ,3]) 
    dev.off()
}

{
    png(sprintf("%s/summaries/state_eff.png", dir))
    par(mfrow=c(1,3))
    state_eff = rstan::extract(fit, "state_eff")[[1]]
    state_eff = apply(state_eff,c(2,3), median)
    hist(state_eff[ ,1])
    hist(state_eff[ ,2])
    hist(state_eff[ ,3]) 
    dev.off()
}

{
    png(sprintf("%s/summaries/scale_state_eff.png", dir))
    par(mfrow=c(1,3))
    hist(rstan::extract(fit, "scale_state_eff[1]")[[1]])
    hist(rstan::extract(fit, "scale_state_eff[2]")[[1]])
    hist(rstan::extract(fit, "scale_state_eff[3]")[[1]])
    dev.off()
}


{
    png(sprintf("%s/summaries/ar_scale.png", dir))
    hist(rstan::extract(fit, "ar_scale")[[1]])
    dev.off()
}


{
    png(sprintf("%s/summaries/timing.png", dir))
    par(mfrow=c(1, 2))
    hist(rstan::extract(fit, "beta_covars_post")[[1]][ , 1, 1])
    hist(rstan::extract(fit, "beta_covars_post")[[1]][ , 1, 2])
    dev.off()
}


{
    png(sprintf("%s/summaries/overdisp.png", dir))
    hist(rstan::extract(fit, "overdisp")[[1]])
    dev.off()
}

if (model_data$bent_cable) {
    {
        png(sprintf("%s/summaries/duration.png", dir))
        hist(6 * rstan::extract(fit, "duration_unc")[[1]])
        dev.off()
    }

    {
        png(sprintf("%s/summaries/lag.png", dir))
        hist(11 + 6 * rstan::extract(fit, "lag_unc")[[1]])
        dev.off()
    }
}
