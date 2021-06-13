library(tidyverse)
library(rstan)
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

fit = read_rds(sprintf("%s/fit.rds", dir))
model_data = read_rds(sprintf("%s/model_data.rds", dir))

ysamples = posterior_predict(
    fit,
    model_data,
    bent_cable=model_data$bent_cable,
    temporal=model_data$temporal,
    spatial=model_data$spatial
)$y_samples
county_lp_var = exp(posterior_predict(
    fit,
    model_data,
    bent_cable=model_data$bent_cable,
    temporal=FALSE,
    spatial=model_data$spatial
)$log_yhat)
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

for (i in 1:length(fiplist)) {
    f1 = fiplist[i]
    ix = which(county_train$fips == f1)
    if (length(ix) == 0)
        next
    pop = 1
    yi = county_train$y[ix] / pop
    predi = ysamples[, ix] / pop
    predi2 = county_lp_var[ ,ix] / pop
    yhati = apply(predi, 2, median)
    yhati2 = apply(predi2, 2, median)
    yhati_95 = apply(predi, 2, quantile, .95)
    yhati_05 = apply(predi, 2, quantile, .05)
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
