library(tidyverse)
library(rstan)
library(abind)
library(reticulate)
library(igraph)

np = import("numpy")
scipy_stats = import("scipy.stats")
nbinom = scipy_stats$nbinom


stan_input_data = function(
  county_train,
  type=c("stayhome", "decrease"),
  lag=12,
  order=2,
  use_mask=FALSE,
  edges=NULL,
  ar_scale=0.25
) {
  county_train = county_train %>%
    mutate(fips_f = as.factor(fips)) %>% 
    mutate(fips_id = as.integer(fips_f)) %>% 
    mutate(state_f = as.factor(state)) %>% 
    mutate(state_id = as.integer(state_f)) %>% 
    arrange(fips_id, days_since_thresh)
  
  fips_ids = levels(county_train$fips_f)
  states_ids = levels(county_train$state_f)
  
  tscale = 100.0
  t1 = county_train$days_since_thresh / tscale
  t1_2 = t1^2
  tpoly_pre = cbind(1, t1)
  if (order > 1)
    for (j in 2:order)
      tpoly_pre = cbind(tpoly_pre, t1^j)  
  
  # this will take time but we want to compute a temporal index
  t = county_train$days_since_thresh
  t = 1 + t - min(t) # so it starts from 1
  fips = county_train$fips
  fips_t_index = setNames(1:nrow(county_train), nm=paste(fips, t, sep="_"))
  fips_tm1 = paste(fips, t - 1, sep="_")
  fips_tm1[t > 1]
  
  k = 1
  k_ = 2
  ar_edges1 = c()
  ar_edges2 = c()
  ar_starts = 1
  for (i in 2:length(t)) {
    if (fips[i - 1] == fips[i]) {
      ar_edges1[k] = i
      ar_edges2[k] = i - 1
      k = k + 1
    } else{
      ar_starts[k_] = i
      k_ = k_ + 1
    }
  }
  
  tm1_pointer = fips_t_index[fips_tm1]
  tm1_pointer[is.na(tm1_pointer)] = 0
  tm1_pointer = tm1_pointer
  
  fips_id = county_train$fips_id
  county_brks = c(0, which(fips_id[-1] != fips_id[-length(fips_id)]), length(fips_id))
  
  
  if (type[1] == "stayhome") {
    days_since_intrv = county_train$days_since_intrv_stayhome
    days_btwn = county_train$days_btwn_stayhome_thresh
  } else if (type[1] == "decrease") {
    days_since_intrv = county_train$days_since_intrv_decrease
    days_btwn = county_train$days_btwn_decrease_thresh
  }
  days_since_intrv[is.na(days_since_intrv)] = -1e6

  t2 = pmax((days_since_intrv - lag) / tscale, 0)
  tpoly_post = matrix(t2, ncol=1)
  if (order > 1)
    for (j in 2:order)
      tpoly_post = cbind(tpoly_post, t2^j)
  
  fips = levels(county_train$fips_f)
  y = county_train$y
  offset = log(county_train$pop) - mean(log(county_train$pop)) - mean(log(1 + y))
  X_pre = as.matrix(select(county_train, college, age_65_plus, black, hispanic))
  X_post = matrix(days_btwn, ncol=1) / tscale
  X_post[is.na(X_post)] = 0.0
  
  mu_X_pre = apply(X_pre, 2, mean)
  sd_X_pre = apply(X_pre, 2, sd)
  for (j in 1:ncol(X_pre))
    X_pre[ ,j] = (X_pre[ ,j] - mu_X_pre[j]) / sd_X_pre[j]
  mu = mu_X_pre
  sig = sd_X_pre
  varname = names(mu_X_pre)
  normalizing_stats = tibble(variable=varname, mean=mu, sd=sig)
  
  # mu_X_post = apply(X_post, 2, mean)
  # sd_X_post = apply(X_post, 2, sd)
  # for (j in 1:ncol(X_post))
  #   X_post[ ,j] = (X_post[ ,j] - mu_X_post[j]) / sd_X_post[j]
  # 
  # mu = c(mu_X_pre, mu_X_post)
  # sig = c(sd_X_pre, sd_X_post)
  # varname = c(names(mu_X_pre), names(mu_X_post))
  # normalizing_stats = tibble(variable=varname, mean=mu, sd=sig)
  # write_csv(normalizing_stats, "normalizing_stats.csv")
  
  time_id = as.integer(county_train$days_since_thresh)
  time_id = 1 + time_id - min(time_id)
  Tmax = max(time_id)
  nchs_id = as.integer(county_train$nchs)
  time_nchs_id = time_id + Tmax * (nchs_id - 1)

  if (!("mask" %in% names(county_train)))
    county_train$mask = rep(1, nrow(county_train))
  
  output = list(
    N = nrow(county_train),
    D_pre = 4,
    D_post = 1,
    fips_ids=fips_ids,
    M = length(fips),
    N_states = length(states_ids),
    X_pre = X_pre,
    X_post = X_post,
    offset = offset,
    tpoly_pre = tpoly_pre,
    days_since_intrv=days_since_intrv,
    y = y,
    tpoly_post=tpoly_post,
    time_id = time_id,
    time_nchs_id = time_nchs_id,
    Tmax = Tmax,
    nchs_id = nchs_id,
    mask = county_train$mask,
    use_mask=use_mask,
    order=order,
    county_id = county_train$fips_id,
    state_id = county_train$state_id,
    normalizing_stats=normalizing_stats,
    tm1_index=tm1_pointer,
    county_brks=county_brks,
    county_lens=table(fips_id),
    ar_edges1=ar_edges1,
    ar_edges2=ar_edges2,
    ar_starts=ar_starts,
    lag=lag,
    df=county_train,
    type=type[1],
    order=order,
    edges=edges,
    ar_scale=ar_scale
  )

  if (!is.null(edges)) {
    # add graph data
    fips = output$fips_ids
    fips2id = setNames(1:length(output$fips_ids), fips)
    dist_scale = 50

    # now the titanic
    edges_ = edges %>%
      mutate(
        wts = ifelse(
          isnbr,
          1.0,
          pmin(1.0, 0.5 * (dist_scale / dist)^2)
        )
      )
    g = igraph::graph_from_data_frame(
      select(edges_, src_lab, tgt_lab),
      vertices=fips,
      directed=FALSE
    )
    output$N_edges = nrow(edges_)
    output$edge_weights = edges_$wts
    comps = components(g)
    output$csizes = comps$csize
    output$csorted = fips2id[names(comps$membership)]
    output$node1 = fips2id[edges_$src_lab]
    output$node2 = fips2id[edges_$tgt_lab]
    output$cmemb = county_train %>% 
      left_join(
        tibble(fips=names(comps$membership), cmemb=comps$membership)
      ) %>% 
      pull(cmemb)
    output$N_comps = comps$no
    output$cbrks = c(0, cumsum(comps$csize))
    output$n_nbrs = map_int(adjacent_vertices(g, fips), length)
    nbrs_eff = edges_ %>% 
      mutate(tmp=src_lab, src_lab=tgt_lab, tgt_lab=tmp) %>% 
      bind_rows(edges_) %>% 
      group_by(src_lab) %>% 
      summarize(nbrs_eff = sum(wts), .groups="drop")

    scale_factor = tibble(
      src_lab = names(comps$membership),
      membership = comps$membership
    ) %>%
      left_join(nbrs_eff, by="src_lab") %>% 
      group_by(membership) %>% 
      summarize(scale_factor =  0.7 * sqrt(mean(nbrs_eff)), .groups="drop")
    output$inv_scaling_factor = 1.0 / scale_factor$scale_factor
    scale_factor$scale_factor[is.na(scale_factor$scale_factor)] = 0
    output$inv_scaling_factor[is.na(output$inv_scaling_factor)] = 0
    output$scaling_factor = scale_factor$scale_factor
    output$nbrs_eff = nbrs_eff$nbrs_eff
    output
  }

  output
}

# posterior predicts outputs the predicted deaths
posterior_predict = function (
  fit,
  model_data,
  states=TRUE,  # if uses tate random effects
  spatial=FALSE,
  rand_lag=FALSE,
  temporal=FALSE,
  cable_bent=FALSE,
  shift_timing=0
) {
  type = model_data$type
  order = model_data$order
  lag = model_data$lag

  parnames = c(
    "nchs_pre",
    "nchs_post",
    "beta_covars_pre",
    "beta_covars_post",
    "baseline_pre",
    "baseline_post",
    "overdisp",
    "rand_eff",
    "scale_rand_eff",
    "log_rate_pre_interv"
  )
  if (spatial)
    parnames = c(parnames, "spatial_eff")
  if (rand_lag)
    parnames = c(parnames, "lag")
  if (cable_bent)
    parnames = c(parnames, c("lag", "duration"))
  if (states)
    parnames = c(parnames, "state_eff")
  if (temporal)
    parnames = c(parnames, "time_term")

  new_df = model_data$df %>% 
    mutate(
      days_since_intrv_stayhome = days_since_intrv_stayhome - shift_timing,
      days_btwn_stayhome_thresh = days_btwn_stayhome_thresh + shift_timing,
      days_since_intrv_decrease = days_since_intrv_decrease - shift_timing,
      days_btwn_decrease_thresh = days_btwn_decrease_thresh + shift_timing,
    )
  # instead of calling input data gain
  # we could just recompute tpoly post down
  # as it is done for rand and cable bent already
  new_data = stan_input_data(new_df, type, order=order, lag=lag, edges=model_data$edges)

  pars = rstan::extract(fit, pars=parnames)
  N = length(new_data$y)
  nsamples = nrow(pars$nchs_pre)
  
  # check compatible size of new data and previous
  tpoly_pre = np$expand_dims(new_data$tpoly_pre, 0L)
  rand_eff_unrolled = np$array(pars$rand_eff[ ,new_data$county_id, ])
  
  if (states) {
    state_eff_unrolled = np$array(pars$state_eff[ ,new_data$state_id, ])
    rand_eff_unrolled = rand_eff_unrolled + state_eff_unrolled
  }
  
    
  rand_eff_term = np$sum(
    np$multiply(tpoly_pre, rand_eff_unrolled), -1L
  )
  
  # eval compatibility of pre_term
  pre_term_prev = pars$log_rate_pre_interv
  size_prev = ncol(pre_term_prev)
  size_new = nrow(new_df)

  eval_pre = shift_timing != 0

  if (eval_pre) {
    X_pre = new_data$X_pre
    X_pre = np$expand_dims(X_pre, 0L)
    X_pre = np$expand_dims(X_pre, 3L)    
    beta_covars_pre = np$expand_dims(np$array(pars$beta_covars_pre), 1L)
    covar_baseline_pre = np$sum(np$multiply(X_pre, beta_covars_pre), -2L)
    nchs_pre_unrolled = np$array(pars$nchs_pre[ ,new_data$nchs_id, ])
    baseline_pre = np$expand_dims(pars$baseline_pre, 1L)
    pre_term = np$add(np$add(baseline_pre, nchs_pre_unrolled), covar_baseline_pre)
    pre_term = np$sum(np$multiply(pre_term, tpoly_pre), -1L)
    offset_ = t(matrix(rep(new_data$offset, times=nsamples), ncol=nsamples))
    pre_term = pre_term + rand_eff_term + offset_
    if (temporal)
      pre_term = pre_term + np$array(pars$time_term)
  } else {
    pre_term = np$array(pars$log_rate_pre_interv)
  }
  
  if (temporal)
    pre_term = pre_term + np$array(pars$time_term)

  # recompute tpoly_post from posterior for post to do counterfactual
  if (rand_lag) {
    lag = np$expand_dims(pars$lag, 1L)
    days_since_intrv_ = np$expand_dims(new_data$days_since_intrv, 0L)
    days_since_intrv_ = as.array(np$add(- lag, days_since_intrv_))
    days_since_intrv_ = 0.01 * pmax(days_since_intrv_, 0)
    tpoly_post = array(days_since_intrv_, dim=c(nsamples, N, 1))
    if (order > 1)
      for (j in 2:order)
        tpoly_post = np$array(abind::abind(tpoly_post, days_since_intrv_^j, along=3))
  } else if (cable_bent) {
    tau = 0.01 * np$expand_dims(pars$lag, 1L)
    gam = 0.01 * np$expand_dims(pars$duration, 1L)
    t_ = 0.01 * np$expand_dims(new_data$days_since_intrv, 0L)
    # transition
    aux1 = np$square(np$maximum(0, np$add(t_, -tau +  gam)))
    aux1 = 0.25 * np$divide(aux1, gam)
    linear1 = aux1 * np$less_equal(t_, np$add(tau, gam))
    aux2 = np$maximum(np$add(t_, -tau), 0)
    linear2 = aux2 * np$greater_equal(t_, np$add(tau, gam))
    linear = linear1 + linear2
    quad = np$square(np$maximum(0, np$add(t_, - tau - gam)))
    tpoly_post = np$stack(list(linear, quad), -1L)
    
  } else {
    tpoly_post = np$expand_dims(new_data$tpoly_post, 0L)
  }

  #
  X_post = new_data$X_post
  X_post = np$expand_dims(X_post, 0L)
  X_post = np$expand_dims(X_post, 3L)
  beta_covars_post = np$expand_dims(np$array(pars$beta_covars_post), 1L)
  covar_baseline_post = np$sum(np$multiply(X_post, beta_covars_post), -2L)
  nchs_post_unrolled = np$array(pars$nchs_post[ ,new_data$nchs_id, , drop=FALSE])
  baseline_post = np$expand_dims(pars$baseline_post, 1L)
  post_term = np$add(
    np$add(baseline_post, nchs_post_unrolled),
    covar_baseline_post
  )
  post_term = np$sum(
    np$multiply(post_term, tpoly_post), -1L
  )
  post_term = as.array(post_term)
  
  log_rate = pre_term + post_term
  rate = exp(log_rate)
  overdisp = matrix(pars$overdisp, nrow=nsamples, ncol=N)
  var = rate + rate ^ 2 / overdisp
  p = pmax((var - rate) / var, 1e-6)
  nbinom_samples = as.array(nbinom$rvs(overdisp, 1 - p))
  #
  list(
    # yhat=rate,
    log_yhat=log_rate,
    # log_yhat_no_rand_eff=log_rate - rand_eff_term,
    # yhat_no_rand_eff=exp(log_rate - rand_eff_term),
    rand_eff=rand_eff_term,
    pre_term=pre_term,
    post_term=post_term,
    y_samples=nbinom_samples
  )
}



# my_posterior_predict_rand_lag = function (
#   fit,
#   new_df,
#   type=c("stayhome", "decrease"),
#   eval_pre = TRUE
# ) {
#   new_data = stan_input_data(new_df, type, lag=0)
#   parnames = c(
#     "nchs_pre", "nchs_post", "beta_covars_pre",
#     "beta_covars_post", "beta_covars_post",
#     "baseline_pre", "baseline_post",
#     "overdisp", "rand_eff",
#     "Omega_rand_eff", "scale_rand_eff",
#     "log_rate_pre_interv", "lag"
#   )
#   pars = rstan::extract(fit, pars=parnames)
#   N = length(new_data$y)
#   nsamples = nrow(pars$nchs_pre)
  
#   # check compatible size of new data and previous
#   tpoly_pre = np$expand_dims(new_data$tpoly_pre, 0L)
#   rand_eff_unrolled = np$array(pars$rand_eff[ ,new_data$county_id, ])
#   rand_eff_term = np$sum(
#     np$multiply(tpoly_pre, rand_eff_unrolled), -1L
#   )
  
#   # eval compatibility of pre_term
#   pre_term_prev = pars$log_rate_pre_interv
#   size_prev = ncol(pre_term_prev)
#   size_new = nrow(new_df)
#   re_eval = eval_pre || (size_prev != size_new)
#   if (re_eval) {
#     X_pre = new_data$X_pre
#     X_pre = np$expand_dims(X_pre, 0L)
#     X_pre = np$expand_dims(X_pre, 3L)    
#     beta_covars_pre = np$expand_dims(np$array(pars$beta_covars_pre), 1L)
#     covar_baseline_pre = np$sum(np$multiply(X_pre, beta_covars_pre), -2L)
#     nchs_pre_unrolled = np$array(pars$nchs_pre[ ,new_data$nchs_id, ])
#     baseline_pre = np$expand_dims(pars$baseline_pre, 1L)
#     pre_term = np$add(np$add(baseline_pre, nchs_pre_unrolled), covar_baseline_pre)
#     pre_term = np$sum(np$multiply(pre_term, tpoly_pre), -1L)
#   } else {
#     pre_term = np$array(pars$log_rate_pre_interv)
#   }

#   #
#   lag = np$expand_dims(pars$lag, 1L)
#   days_since_intrv_ = np$expand_dims(new_data$days_since_intrv, 0L)
#   days_since_intrv_ = as.array(np$add(- lag, days_since_intrv_))
#   days_since_intrv_ = 0.01 * pmax(days_since_intrv_, 0) 
#   tpoly_post = abind::abind(days_since_intrv_, days_since_intrv_^2, along=3)
#   #
#   X_post = new_data$X_post
#   X_post = np$expand_dims(X_post, 0L)
#   X_post = np$expand_dims(X_post, 3L)
#   beta_covars_post = np$expand_dims(np$array(pars$beta_covars_post), 1L)
#   covar_baseline_post = np$sum(np$multiply(X_post, beta_covars_post), -2L)
#   nchs_post_unrolled = np$array(pars$nchs_post[ ,new_data$nchs_id, ])
#   baseline_post = np$expand_dims(pars$baseline_post, 1L)
#   post_term = np$add(
#     np$add(baseline_post, nchs_post_unrolled),
#     covar_baseline_post
#   )
#   post_term = np$sum(
#     np$multiply(post_term, tpoly_post), -1L
#   )
#   post_term = as.array(post_term)
  
#   log_rate = pre_term + post_term
#   rate = exp(log_rate)
#   overdisp = matrix(pars$overdisp, nrow=nsamples, ncol=N)
#   var = rate + rate ^ 2 / overdisp
#   p = pmax((var - rate) / var, 1e-6)
#   nbinom_samples = as.array(nbinom$rvs(overdisp, 1 - p))
#   #
#   list(
#     yhat=rate,
#     tpoly_post=tpoly_post,
#     log_yhat=log_rate,
#     log_yhat_no_rand_eff=log_rate - rand_eff_term,
#     yhat_no_rand_eff=exp(log_rate - rand_eff_term),
#     rand_eff_term=rand_eff_term,
#     pre_term=pre_term,
#     post_term=post_term,
#     y_samples=nbinom_samples
#   )
# }



# my_posterior_predict_cable_bent = function (
#   fit,
#   new_df,
#   type=c("stayhome", "decrease"),
#   eval_pre = TRUE
# ) {
#   new_data = stan_input_data(new_df, type, lag=0)
#   parnames = c(
#     "nchs_pre", "nchs_post", "beta_covars_pre",
#     "beta_covars_post", "beta_covars_post",
#     "baseline_pre", "baseline_post",
#     "overdisp", "rand_eff",
#     "Omega_rand_eff", "scale_rand_eff",
#     "log_rate_pre_interv", "lag", "duration"
#   )
#   pars = rstan::extract(fit, pars=parnames)
#   N = length(new_data$y)
#   nsamples = nrow(pars$nchs_pre)
  
#   # check compatible size of new data and previous
#   tpoly_pre = np$expand_dims(new_data$tpoly_pre, 0L)
#   rand_eff_unrolled = np$array(pars$rand_eff[ ,new_data$county_id, ])
#   rand_eff_term = np$sum(
#     np$multiply(tpoly_pre, rand_eff_unrolled), -1L
#   )
  
#   # eval compatibility of pre_term
#   pre_term_prev = pars$log_rate_pre_interv
#   size_prev = ncol(pre_term_prev)
#   size_new = nrow(new_df)
#   re_eval = eval_pre || (size_prev != size_new)
#   if (re_eval) {
#     X_pre = new_data$X_pre
#     X_pre = np$expand_dims(X_pre, 0L)
#     X_pre = np$expand_dims(X_pre, 3L)    
#     beta_covars_pre = np$expand_dims(np$array(pars$beta_covars_pre), 1L)
#     covar_baseline_pre = np$sum(np$multiply(X_pre, beta_covars_pre), -2L)
#     nchs_pre_unrolled = np$array(pars$nchs_pre[ ,new_data$nchs_id, ])
#     baseline_pre = np$expand_dims(pars$baseline_pre, 1L)
#     pre_term = np$add(np$add(baseline_pre, nchs_pre_unrolled), covar_baseline_pre)
#     pre_term = np$sum(np$multiply(pre_term, tpoly_pre), -1L)
#     offset_ = t(matrix(rep(new_data$offset, times=nsamples), ncol=nsamples))
#     pre_term = pre_term + rand_eff_term + offset_
#   } else {
#     pre_term = np$array(pars$log_rate_pre_interv)
#   }
  
#   #
#   lag = 0.01 * np$expand_dims(pars$lag, 1L)
#   t_ = 0.01 * np$expand_dims(new_data$days_since_intrv, 0L)
#   gam = 0.01 * np$expand_dims(pars$duration, 1L)

#   transition = as.array(
#     0.25 * np$divide(np$square(np$add(t_, -lag + gam)), gam)
#     * np$less_equal(t_, lag + gam)
#     * np$greater_equal(t_, lag - gam)
#   )
#   post1 = as.array(
#     np$subtract(t_, lag) * np$greater_equal(t_, lag + gam)
#   )
#   post2 = as.array(
#     np$square(np$maximum(np$subtract(t_, lag + gam), 0))
#   )
#   tpoly_post = abind::abind(transition + post1, post2, along=3)
#   #
#   X_post = new_data$X_post
#   X_post = np$expand_dims(X_post, 0L)
#   X_post = np$expand_dims(X_post, 3L)
#   beta_covars_post = np$expand_dims(np$array(pars$beta_covars_post), 1L)
#   covar_baseline_post = np$sum(np$multiply(X_post, beta_covars_post), -2L)
#   nchs_post_unrolled = np$array(pars$nchs_post[ ,new_data$nchs_id, ])
#   baseline_post = np$expand_dims(pars$baseline_post, 1L)
#   post_term = np$add(
#     np$add(baseline_post, nchs_post_unrolled),
#     covar_baseline_post
#   )
#   post_term = np$sum(
#     np$multiply(post_term, tpoly_post), -1L
#   )
#   post_term = as.array(post_term)
  
#   log_rate = pre_term + post_term
#   rate = exp(log_rate)
#   overdisp = matrix(pars$overdisp, nrow=nsamples, ncol=N)
#   var = rate + rate ^ 2 / overdisp
#   p = pmax((var - rate) / var, 1e-6)
#   nbinom_samples = as.array(nbinom$rvs(overdisp, 1 - p))
#   #
#   list(
#     yhat=rate,
#     tpoly_post=tpoly_post,
#     log_yhat=log_rate,
#     log_yhat_no_rand_eff=log_rate - rand_eff_term,
#     yhat_no_rand_eff=exp(log_rate - rand_eff_term),
#     rand_eff_term=rand_eff_term,
#     pre_term=pre_term,
#     post_term=post_term,
#     y_samples=nbinom_samples
#   )
# }
