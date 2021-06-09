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
  lag=14,
  order=2,
  use_mask=FALSE,
  edges=NULL,
  ar_scale=0.25, 
  old_model_data=NULL,
  all_post=FALSE,  # it's a legacy shortcut, should be removed
  pre_vars = c("college", "age_65_plus", "black", "hispanic", "popdensity"),
  post_vars = c("college", "age_65_plus", "black", "hispanic", "popdensity"),
  pre_inter_vars = c("college", "age_65_plus", "black", "hispanic", "popdensity"),
  post_inter_vars = c(),
  use_post_inter = FALSE,
  use_pre_inter = TRUE,
  spatial_scale_fixed = 1.0,
  autocor=0.7,
  k_nbrs=2,
  ar_tight_prior_scale=FALSE
) {
  
  if (is.null(old_model_data)) {
    fips_ids = unique(county_train$fips)
    states_ids = unique(county_train$state)
  } else{
    fips_ids = old_model_data$fips_ids
    states_ids = old_model_data$states_ids
    lag = old_model_data$lag
    use_post_inter = old_model_data$lag
    use_pre_inter = old_model_data$lag
    pre_inter_vars = old_model_data$pre_inter_vars
    post_inter_vars = old_model_data$post_inter_vars
    pre_vars = old_model_data$pre_vars
    post_vars = old_model_data$post_vars
    spatial_scale_fixed = old_model_data$spatial_scale_fixed
    ar_scale = old_model_data$ar_scale
    autocor = old_model_data$autocor
  }
  
  county_train = county_train %>%
    mutate(fips_f = factor(fips, levels=fips_ids)) %>% 
    mutate(fips_id = as.integer(fips_f)) %>% 
    mutate(state_f = factor(state, levels=states_ids)) %>% 
    mutate(state_id = as.integer(state_f)) %>% 
    arrange(fips_id, days_since_thresh)
  
  
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
  
  if (is.null(old_model_data)) {
    y = county_train$y
    offset_baseline = - mean(log(1 + y))
  } else{
    y = old_model_data$y
    offset_baseline = old_model_data$offset_baseline
  }
  
  X_pre = as.matrix(county_train[ ,pre_vars])
  X_pre_inter = as.matrix(county_train[ ,pre_inter_vars])

  X_post_inter = 0.01 * matrix(county_train$days_since_thresh, ncol=1)
  if (length(post_inter_vars) > 1)
    X_post_inter = cbind(X_post_inter, as.matrix(county_train[ ,post_inter_vars]))
  
  offset = log(county_train$pop) - mean(log(county_train$pop)) + offset_baseline
  
  if (is.null(old_model_data)) {
    mu_X_pre = apply(X_pre, 2, mean)
    sd_X_pre = apply(X_pre, 2, sd)
  } else{
    mu_X_pre = old_model_data$normalizing_stats$mean
    sd_X_pre = old_model_data$normalizing_stats$sd
  }
  # print(dim(X_pre))
  for (j in 1:ncol(X_pre))
    X_pre[ ,j] = (X_pre[ ,j] - mu_X_pre[j]) / sd_X_pre[j]
  mu = mu_X_pre
  sig = sd_X_pre
  varname = names(mu_X_pre)
  normalizing_stats = tibble(variable=varname, mean=mu, sd=sig)
  
  X_post = matrix(days_btwn, ncol=1) / tscale
  X_post[is.na(X_post)] = 0.0
  X_post = cbind(X_post, as.matrix(county_train[ ,post_vars]))
  
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
  
  X_post_inter = X_post
  X_pre_inter = X_pre
  
  output = list(
    N = nrow(county_train),
    D_pre = ncol(X_pre),
    D_pre_inter = ncol(X_pre_inter),
    D_post = ncol(X_post),
    D_post_inter = ncol(X_post_inter),
    fips_ids=fips_ids,
    states_ids=states_ids,
    M = length(fips),
    N_states = length(states_ids),
    X_pre = X_pre,
    X_pre_inter = X_pre_inter,
    X_post = X_post,
    X_post_inter = X_post_inter,
    offset = offset,
    offset_baseline=offset_baseline,
    tpoly_pre = tpoly_pre,
    days_since_intrv=days_since_intrv,
    y = y,
    tpoly_post=tpoly_post,
    time_id = time_id,
    time_nchs_id = time_nchs_id,
    Tmax = Tmax,
    nchs_id = nchs_id,
    mask = county_train$mask,
    mask_obs = which(as.logical(county_train$mask)),
    mask_miss = which(!as.logical(county_train$mask)),
    N_miss = sum(!as.logical(county_train$mask)),
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
    ar_tight_prior_scale=as.integer(ar_tight_prior_scale),
    lag=lag,
    df=county_train,
    type=type[1],
    order=order,
    edges=edges,
    ar_scale=ar_scale,
    use_post_inter=as.numeric(use_post_inter),
    use_pre_inter=as.numeric(use_pre_inter),
    pre_vars=pre_vars,
    post_vars=post_vars,
    pre_inter_vars=pre_inter_vars,
    post_inter_vars=post_inter_vars,
    spatial_scale_fixed=spatial_scale_fixed,
    autocor=autocor
  )
  
  if (!is.null(edges)) {
    # add graph data
    fips = output$fips_ids
    fips2id = setNames(1:length(output$fips_ids), fips)
    dist_scale = 50

    # make nearest neighbors graph
    edges_ = edges %>% 
      ungroup() %>% 
      mutate(tmp=src_lab, src_lab=tgt_lab, tgt_lab=tmp, tmp=src, src=tgt, tgt=tmp) %>% 
      bind_rows(edges) %>%
      arrange(src, dist) %>% 
      group_by(src) %>% 
      slice(1:k_nbrs) %>%  # keep 2 nearest neighbors
      ungroup() %>%
      bind_rows(filter(edges, isnbr)) %>%
      mutate(src_lab_ = pmin(src_lab, tgt_lab), tgt_lab_ = pmax(src_lab, tgt_lab)) %>% 
      distinct(src_lab_, tgt_lab_, .keep_all = TRUE) %>% 
      select(-tmp)

    # # now the titanic
    edges_ = edges_ %>%
      mutate(
        wts = ifelse(
          isnbr,
          1.0,
          1.0 # pmin(1.0, (dist_scale / dist))
        )
      )
    output$edges_post = edges_
    g = igraph::graph_from_data_frame(
      select(edges_, src_lab, tgt_lab),
      vertices=unique(fips),
      directed=FALSE
    )
    output$n_nbrs = map_int(adjacent_vertices(g, fips), length)
    output$N_edges = nrow(edges_)
    output$edge_weights = edges_$wts
    comps = components(g)
    output$csizes = comps$csize
    output$csorted = fips2id[names(comps$membership)]
    atleast2 = output$n_nbrs > 2
    output$node1 = fips2id[edges_$src_lab]
    output$node2 = fips2id[edges_$tgt_lab]
    # 
    output$cmemb = county_train %>%
      left_join(
        tibble(fips=names(comps$membership), cmemb=comps$membership)
      ) %>%
      pull(cmemb)
    output$N_comps = comps$no
    output$cbrks = c(0, cumsum(comps$csize))

    nbrs_eff = edges_ %>%
      mutate(tmp=src_lab, src_lab=tgt_lab, tgt_lab=tmp) %>%
      bind_rows(edges_) %>%
      group_by(src_lab) %>%
      summarize(nbrs_eff = sum(wts), nbrs_d = length(wts), .groups="drop")

    scale_factor = tibble(
      src_lab = names(comps$membership),
      membership = comps$membership
    ) %>%
      left_join(nbrs_eff, by="src_lab") %>% 
      group_by(membership) %>% 
      summarize(mean_nbrs_d = mean(nbrs_d), mean_nbrs_eff = mean(nbrs_eff), scale_factor =  0.7 * sqrt(mean(nbrs_eff)), .groups="drop")
    
    
    scale_factor2 = tibble(
      src_lab = names(comps$membership),
      membership = comps$membership
    )  %>%
      left_join(nbrs_eff, by="src_lab") %>% 
      group_by(membership) %>% 
      mutate(mean_nbrs = rep(mean(nbrs_d), n()), mean_nbrs_eff = rep(mean(nbrs_eff), n()), csize = n()) %>% 
      ungroup()
    mean_nbrs = setNames(scale_factor2$mean_nbrs, scale_factor2$src_lab)
    mean_nbrs_eff = setNames(scale_factor2$mean_nbrs_eff, scale_factor2$src_lab)
    
    output$cmeannbrs_eff = scale_factor$mean_nbrs_eff[output$cmemb]
    output$cmeannbrs_eff[is.na(output$cmeannbrs_eff)] = 0
    output$cmeannbrs_d = scale_factor$mean_nbrs_d[output$cmemb]
    output$cmeannbrs_d[is.na(output$cmeannbrs_d)] = 0
    output$inv_scaling_factor = 1.0 / scale_factor$scale_factor
    scale_factor$scale_factor[is.na(scale_factor$scale_factor)] = 0
    output$inv_scaling_factor[is.na(output$inv_scaling_factor)] = 0
    output$scaling_factor = scale_factor$scale_factor
    output$nbrs_eff = nbrs_eff$nbrs_eff
    output$bym_scaled_edge_weights = output$edge_weights / (0.7 * sqrt(mean_nbrs_eff[output$node1]))
    output
  }
  
  output
}

# posterior predicts outputs the predicted deathss
posterior_predict = function (
  fit,
  model_data,
  new_df=NULL,
  rand_eff=TRUE,
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
  pre_inter = as.logical(model_data$use_pre_inter)
  post_inter = as.logical(model_data$use_post_inter)
  
  parnames = c(
    "nchs_pre",
    "nchs_post",
    "beta_covars_pre",
    "beta_covars_post",
    "baseline_pre",
    "baseline_post",
    "overdisp",
    "rand_eff",
    "scale_rand_eff"
  )
  
  if (spatial)
    parnames = c(parnames, "spatial_eff")
  if (rand_lag)
    parnames = c(parnames, "lag_unc")
  if (cable_bent)
    parnames = c(parnames, c("lag_unc", "duration_unc"))
  if (states)
    parnames = c(parnames, "state_eff")
  if (temporal)
    parnames = c(parnames, "time_term")
  if (pre_inter)
    parnames = c(parnames, "beta_covars_pre_inter")
  if (post_inter)
    parnames = c(parnames, "beta_covars_post_inter")
  
  if (is.null(new_df))
    new_df = model_data$df
  
  new_df = new_df %>% 
    mutate(
      days_since_intrv_stayhome = days_since_intrv_stayhome - shift_timing,
      days_btwn_stayhome_thresh = days_btwn_stayhome_thresh + shift_timing,
      days_since_intrv_decrease = days_since_intrv_decrease - shift_timing,
      days_btwn_decrease_thresh = days_btwn_decrease_thresh + shift_timing,
    )
  # instead of calling input data gain
  # we could just recompute tpoly post down
  # as it is done for rand and cable bent already
  new_data = stan_input_data(
    new_df,
    type,
    edges=model_data$edges,
    old_model_data = model_data
  )
  
  pars = rstan::extract(fit, pars=parnames)
  N = nrow(new_df)
  nsamples = nrow(pars$nchs_pre)
  
  # check compatible size of new data and previous
  tpoly_pre = np$expand_dims(new_data$tpoly_pre, 0L)
  rand_eff_unrolled = np$array(pars$rand_eff[ ,new_data$county_id, ])
  
  rand_eff_term = np$sum(
    np$multiply(tpoly_pre, rand_eff_unrolled), -1L
  )
  
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
  
  pre_term = pre_term + offset_
  
  if (temporal) {
    time_eff = pars$time_term
    M = new_data$M
    brks = new_data$county_brks
    for (j in 1:M) {
      start = brks[j] + 1
      finish = brks[j + 1]
      time_eff[ ,start] = - apply(time_eff[ ,(start + 1):finish], 1, sum)
    }
    pre_term = pre_term + time_eff
  }
  
  if (states) {
    state_eff_unrolled = np$array(pars$state_eff[ ,new_data$state_id, ])
    state_eff_term = np$sum(
      np$multiply(tpoly_pre, state_eff_unrolled), -1L
    )
    pre_term = pre_term + state_eff_term
  }
  
  if (spatial) {
    tpoly_pre = np$expand_dims(new_data$tpoly_pre, 0L)
    spatial_eff_unrolled = np$array(pars$spatial_eff[ ,new_data$county_id, ])
    spatial_eff_term = np$sum(
      np$multiply(tpoly_pre, spatial_eff_unrolled), -1L
    )
    pre_term = pre_term + spatial_eff_term
  }
  
  
  # recompute tpoly_post from posterior for post to do counterfactual
  if (rand_lag) {
    lag = 11 + 5 * np$expand_dims(pars$lag_unc, 1L)
    days_since_intrv_ = np$expand_dims(new_data$days_since_intrv, 0L)
    days_since_intrv_ = as.array(np$add(- lag, days_since_intrv_))
    days_since_intrv_ = 0.01 * pmax(days_since_intrv_, 0)
    tpoly_post = array(days_since_intrv_, dim=c(nsamples, N, 1))
    if (order > 1)
      for (j in 2:order)
        tpoly_post = np$array(abind::abind(tpoly_post, days_since_intrv_^j, along=3))
  } else if (cable_bent) {
    lag = 11 + 5 * np$expand_dims(pars$lag_unc, 1L)
    duration = pars$duration_unc
    tau = 0.01 * np$expand_dims(lag, 1L)
    gam = 0.01 * np$expand_dims(duration, 1L)
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
  
  if (post_inter) {
    days_since_intrv_ = np$expand_dims(new_data$days_since_intrv, 0L)
    days_since_intrv_ = as.array(np$add(- lag, days_since_intrv_))
    days_since_intrv_ = 0.01 * pmax(days_since_intrv_, 0)
    X_post_inter = new_data$X_post_inter
    X_post_inter = np$expand_dims(X_post_inter, 0L)
    X_post_inter = np$expand_dims(X_post_inter, 3L)
    beta_covars_post_inter =  np$expand_dims(pars$beta_covars_post_inter[, new_data$nchs_id, , drop=FALSE], -1L)
    post_inter_term = np$sum(np$multiply(beta_covars_post_inter, X_post_inter), 2L)
    post_inter_term = post_inter_term[, , 1]
    t1 = np$expand_dims(new_data$tpoly_post[ , 1], 0L)
    post_inter_term = np$multiply(t1, post_inter_term)
    post_term = post_term + as.array(post_inter_term)
  }
  
  if (pre_inter) {
    X_pre_inter = new_data$X_pre_inter
    X_pre_inter = np$expand_dims(X_pre_inter, 0L)
    X_pre_inter = np$expand_dims(X_pre_inter, 3L)
    D_inter = new_data$D_pre_inter

    beta_covars_pre_inter0 =  np$expand_dims(pars$beta_covars_pre_inter[, new_data$nchs_id, 1:D_inter, drop=FALSE], -1L)
    beta_covars_pre_inter1 =  np$expand_dims(pars$beta_covars_pre_inter[, new_data$nchs_id, (D_inter + 1):(2 * D_inter), drop=FALSE], -1L)
    beta_covars_pre_inter2 =  np$expand_dims(pars$beta_covars_pre_inter[, new_data$nchs_id, (2 * D_inter + 1):(3 * D_inter), drop=FALSE], -1L)
   
     pre_inter_term0 = np$sum(np$multiply(beta_covars_pre_inter0, X_pre_inter), 2L)[ , , 1]
    
    t1 = np$expand_dims(new_data$tpoly_pre[ , 2], 0L)
    pre_inter_term1 = np$sum(np$multiply(beta_covars_pre_inter1, X_pre_inter), 2L)[ , , 1]
    pre_inter_term1 = np$multiply(t1, pre_inter_term1)
    
    t2 = np$expand_dims(new_data$tpoly_pre[ ,3], 0L)
    pre_inter_term2 = np$sum(np$multiply(beta_covars_pre_inter2, X_pre_inter), 2L)[ , , 1]
    pre_inter_term2 = np$multiply(t2, pre_inter_term2)
  
    pre_term = pre_term + pre_inter_term0 + pre_inter_term1 + pre_inter_term2
  }
  
  log_rate = pre_term + post_term
  if (rand_eff)
    log_rate = log_rate + rand_eff_term
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
