// The input data is a vector 'y' of length 'N'.
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of obs
  int<lower=0> D_pre; // number of covariates before internvention
  int<lower=0> D_post; // number of covariates after internvention
  int<lower=0> D_pre_inter; // number of covariates before internvention
  int<lower=0> D_post_inter; // number of covariates after internvention
  int<lower=0> M; // number of counties
  int<lower=0> y[N]; // deaths
  int<lower=0, upper=N> N_miss;
  int<lower=0, upper=1> mask[N];
  int<lower=1,upper=N> mask_miss[N_miss]; // ignores the y-values for these entries
  int<lower=1,upper=N> mask_obs[N - N_miss]; //
  int<lower=0,upper=1> use_mask; // ignore some data entries (used for cross-validation training)
  real<lower=0, upper=1> use_post_inter;
  real<lower=0, upper=1> use_pre_inter;
  matrix[N, D_pre] X_pre; // covars for pre-trend
  matrix[N, D_pre_inter] X_pre_inter; // covars for pre-trend
  matrix[N, D_post] X_post;  // covars for post-trend, the only one here is days since intervention
  matrix[N, D_post_inter] X_post_inter;  // covars for post-trend, the only one here is days since intervention
  vector[N] offset;  // log of population
  int<lower=0> county_id[N];  // county indicator
  int<lower=0> nchs_id[N];  // nchs indicator
  matrix[N, 3] tpoly_pre; // days since threhsold
  matrix[N, 2] tpoly_post; // days since intervention
  int<lower=0> N_states; // number of states
  int<lower=0> state_id[N];  // county indicator
  //
  // Needed for ICAR prior, it will just penalize edge-wise differneces in the spatial effects
  real<lower=0> spatial_scale_fixed;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  vector<lower=0>[N_edges] edge_weights;
  vector<lower=0>[N_edges] bym_scaled_edge_weights;  // each edge is w_ij / 0.7 sqrt(m) 
                        // where m is the average number of nbrs in the connected
                        // component where the edge is.
  int<lower=0> N_comps;
  // int cmemb[N];   // indicator of connected component membership
  int csizes[N_comps];   // size of each connected component
  // int csorted[M];   // pointer to nodes where they are sorted by component
  int cbrks[N_comps + 1];   // where each component begins and ends in the above pointer
  // for ar(1)
  vector[N] days_since_intrv; // days since/before intervention
  int<lower=1,upper=N> ar_edges1[N - M];
  int<lower=1,upper=N> ar_edges2[N - M];
  int<lower=1,upper=N> ar_starts[M];
  int<lower=0,upper=N> county_brks[M + 1];
  int<lower=0,upper=N> county_lens[M];
  // real spatial_scale;
  // real<lower=0> ar_scale;
  real<lower=0.0, upper=1.0> autocor;
  int<lower=0, upper=1> ar_tight_prior_scale;
}

transformed data{
  int order = 2;  // polynomial order for the trends
  // real acor_mu = 0.95;  //
  // real<lower=0> acor_prec = 2.0 * M;
  // real acor_mu = 0.7;  //
  // real<lower=0> acor_prec = 10.0 * M;
}

parameters {
  matrix<lower=-10.0,upper=10.0>[6, 2] nchs_pre_lin;
  matrix<lower=-10.0,upper=0.0>[6, 1] nchs_pre_quad;
  matrix<lower=-10.0,upper=10.0>[D_pre, order + 1] beta_covars_pre;
  row_vector<lower=-10.0,upper=10.0>[order + 1] baseline_pre;
  matrix<lower=-10.0,upper=10.0>[D_post, order] beta_covars_post;
  matrix<lower=-10.0,upper=10.0>[6, D_post_inter] beta_covars_post_inter;
  matrix<lower=-10.0,upper=10.0>[6, 3 * D_pre_inter] beta_covars_pre_inter;
  row_vector<lower=-10.0,upper=10.0>[order] baseline_post;
  matrix<lower=-10.0,upper=10.0>[6, 1] nchs_post_lin;
  matrix<lower=-10.0,upper=0.0>[6, 1] nchs_post_quad;
  real<lower=0.025, upper=300.0> overdisp;
  //
  corr_matrix[3] Omega_rand_eff;
  corr_matrix[3] Omega_state_eff;
  matrix<lower=-10.0,upper=10.0>[M, 2] rand_eff_lin;
  matrix<lower=-10.0,upper=0.0>[M, 1] rand_eff_quad;
  matrix<lower=-10.0,upper=10.0>[N_states, 2] state_eff_lin;
  matrix<lower=-10.0,upper=0.0>[N_states, 1] state_eff_quad;
  vector<lower=0.001, upper=20.0>[3] scale_rand_eff;
  vector<lower=0.001, upper=20.0>[3] scale_state_eff;
  //
  matrix<lower=-10.0,upper=10.0>[M, 2] spatial_eff_lin;   // scaled spatial effects
  matrix<lower=-10.0,upper=0.0>[M, 1] spatial_eff_quad;   // scaled spatial effects
  vector<lower=0.0, upper=30.0>[3] spatial_scale;

  // for temporal error
  vector<lower=-10.0, upper=10.0>[N] time_term;
  // real<lower=0.0, upper=1.0> autocor;
  real<lower=0> ar_scale;
}

transformed parameters {
  matrix[M, 3] rand_eff = append_col(rand_eff_lin, rand_eff_quad);
  matrix[N_states, 3] state_eff = append_col(state_eff_lin, state_eff_quad);
  matrix[M, 3] spatial_eff = append_col(spatial_eff_lin, spatial_eff_quad);
  matrix[6, 3] nchs_pre = append_col(nchs_pre_lin, nchs_pre_quad);
  matrix[6, 2] nchs_post = append_col(nchs_post_lin, nchs_post_quad);

  vector[N] rand_eff_term = rows_dot_product(
    rand_eff[county_id, :],  // random effects unfolded
    tpoly_pre
  );
  
  vector[N] state_eff_term = rows_dot_product(
    state_eff[state_id, :],  // random effects unfolded
    tpoly_pre
  );

  // spatial efects unrolled
  vector[N] spatial_eff_term = rows_dot_product(
    spatial_eff[county_id, 1:3], tpoly_pre
  );

  vector[N] pre_term = rows_dot_product(
    (
      + X_pre * beta_covars_pre  // interaction with covariates pre-interv
      + rep_matrix(baseline_pre, N)
      + nchs_pre[nchs_id, :]
    ),
    tpoly_pre
  ) + use_pre_inter * (
    rows_dot_product(X_pre_inter, beta_covars_pre_inter[nchs_id, 1:D_pre_inter]) +
    col(tpoly_pre, 2) .* rows_dot_product(X_pre_inter, beta_covars_pre_inter[nchs_id, (D_pre_inter + 1):(2 * D_pre_inter)]) +
    col(tpoly_pre, 3) .* rows_dot_product(X_pre_inter, beta_covars_pre_inter[nchs_id, (2 * D_pre_inter + 1):(3 * D_pre_inter)])
  );
  // if (use_pre_inter == 1)
  //   pre_term = pre_term + pre_inter_term0 + pre_inter_term1;
  vector[N] post_inter_term = col(tpoly_pre, 1) .* rows_dot_product(X_post_inter, beta_covars_post_inter[nchs_id, :]);
  vector[N] post_term = rows_dot_product(
    (
      X_post * beta_covars_post  // interaction with covariates post-interv
      + rep_matrix(baseline_post, N)
      + nchs_post[nchs_id, :]
    ),
    tpoly_post
  ) + use_post_inter * post_inter_term;

  // vector[N] log_rate_pre_interv = offset + rand_eff_term + pre_term + state_eff_term + spatial_eff_term;  # + time_term;
  vector[N] log_rate_pre_interv = offset + rand_eff_term + pre_term + state_eff_term + spatial_eff_term + time_term;
  vector[N] log_rate = log_rate_pre_interv + post_term;
  
  matrix[3, 3] Sigma_state_eff = quad_form_diag(Omega_state_eff, scale_state_eff);
  matrix[3, 3] Sigma_rand_eff = quad_form_diag(Omega_rand_eff, scale_rand_eff);
  // real autocor = autocor_unc;  // for backward compaibility
}

model {
  // parameter priors
  overdisp ~ exponential(1.0);
  // scale_rand_eff ~ normal(0, 1.0 / 17.0);
  scale_rand_eff ~ normal(0, 1.0);
  Omega_rand_eff ~ lkj_corr(2.0);
  scale_state_eff ~ normal(0, 10.0);
  Omega_state_eff ~ lkj_corr(2.0);
  // scale_state_eff ~ normal(0, 1.0);
  // spatial
  // spatial_scale ~ normal(0.0, 2.0 / 17.0);
  spatial_scale ~ normal(0.0, 1.0);
  // spatial_scale ~ normal(0.0, 0.001);

  to_vector(beta_covars_pre) ~ normal(0, 10.0);
  to_vector(beta_covars_post) ~ normal(0, 10.0);
  beta_covars_pre_inter[1, :] ~ normal(0, 0.001);
  beta_covars_post_inter[1, :] ~ normal(0, 0.001);
  to_vector(beta_covars_pre_inter) ~ normal(0, 10.0);
  to_vector(beta_covars_post_inter) ~ normal(0, 10.0);
  nchs_pre[1,:] ~ normal(0, 0.001);
  nchs_post[1,:] ~ normal(0, 0.001);
  to_vector(nchs_pre[2:6, :]) ~ normal(0, 10.0);
  to_vector(nchs_post[2:6, :]) ~ normal(0, 10.0);
  baseline_post ~ normal(0.0, 10.0);
  baseline_pre ~ normal(0.0, 10.0);
  // autocor ~ normal(0.0, )

  // random effect priors (independent)
  for (i in 1:N_states)
    row(state_eff, i) ~ multi_normal(rep_vector(0.0, order + 1), Sigma_state_eff);
  for (i in 1:M)
    row(rand_eff, i) ~ multi_normal(rep_vector(0.0, order + 1), Sigma_rand_eff);
  for (j in 1:(order + 1)) {
    col(state_eff, j) ~ normal(0, 10.0);  // tiny reg
    col(rand_eff, j) ~ normal(0, 10.0);  // tiny reg
  }

  //
  // scale_spatial_eff ~ normal(0, 10.0);
  for (j in 1:3) {
    if (spatial_scale_fixed == 0.0) {
      target += -0.5 * dot_self(bym_scaled_edge_weights .* (spatial_eff[node1, j] - spatial_eff[node2, j])) ./ square(spatial_scale[j]);
    } else {
      target += -0.5 * dot_self(bym_scaled_edge_weights .* (spatial_eff[node1, j] - spatial_eff[node2, j])) ./ square(spatial_scale_fixed);
    }    
    // target += -0.5 * dot_self(edge_weights .* (spatial_eff[node1, j] - spatial_eff[node2, j])) ./ square(spatial_scale);
    for (c in 1:N_comps)
      sum(spatial_eff[(cbrks[c] + 1):cbrks[c], j]) ~ normal(0, 0.001 * csizes[c]);
    col(spatial_eff, j) ~ normal(0, 10.0);  // tiny reg
  }

  // deal with missing obs
  if (ar_tight_prior_scale == 1) {
    ar_scale ~ normal(0.0, 1.0 / sqrt(M));
  } else {
    // ar_scale ~ normal(0.0, 1.0 / sqrt(N));
    // ar_scale ~ normal(0.0, 1.0 / sqrt(M));
    ar_scale ~ normal(0.0, 1.0);
  }

  if (autocor < 1.0) {
    time_term[ar_edges1] ~ normal(autocor * time_term[ar_edges2], ar_scale);
    time_term[ar_starts] ~ normal(0, ar_scale / sqrt(1.0 - square(autocor)));
  } else {
    target += - 0.5 * dot_self(time_term[ar_edges1] - time_term[ar_edges2]) / square(ar_scale);
    time_term[ar_starts] ~ normal(0, ar_scale);
  }

  for (j in 1:M) {
    sum(time_term[(county_brks[j] + 1):(county_brks[j + 1])]) ~ normal(0.0, 0.001 * county_lens[j]);
  }
 
  if (use_mask == 0) {
    y ~ neg_binomial_2(exp(log_rate) + 1e-8, overdisp);
  }
  else {
    y[mask_obs] ~ neg_binomial_2(exp(log_rate[mask_obs]) + 1e-8, overdisp);
    if (autocor < 1.0) {
      time_term[mask_miss] ~ normal(0, 0.1 * ar_scale / sqrt(1.0 - square(autocor))); // this ones aren't observed
    } else {
      time_term[mask_miss] ~ normal(0, 0.1 * ar_scale); // this ones aren't observed
    }
  }


  // if (use_mask == 1)
  //   time_term[mask_miss] ~ normal(0, 0.01); // these ones aren't observed
  // time_term[ar_starts] ~ normal(0, 0.001);  // these ones are fake, never used in likelihood
  // time_term ~ normal(0, 10.0);  // reg  for stability
  // // // for likelihood replace starts with - sum of all
  // for (j in 1:M) {  // for every county
  //   int start = county_brks[j] + 1;
  //   int finish = county_brks[j + 1];
  //   // for t0 ignore the parameter and treat it as the negative sum of all others to guarantee zero-sum
  //   // real sum_all_but_first = sum(time_term[(start + 1):finish]);
  //   if (autocor < 1) {
  //      sum(time_term[(start + 1):finish]) ~ normal(0.0, ar_scale / sqrt(1.0 - square(autocor)));
  //   } else {
  //      sum(time_term[(start + 1):finish]) ~ normal(0.0, ar_scale);      
  //   }
  //   y[start] ~ neg_binomial_2(exp(log_rate[start] - sum(time_term[(start + 1):finish])) + 1e-8, overdisp);
  //   // for the rest do usual ar(1) but make sure tie it with t0 correctly
  //   for (i in (start + 1):finish) {
  //     if (use_mask == 0 || mask[i] == 1)
  //         y[i] ~ neg_binomial_2(exp(log_rate[i] + time_term[i]) + 1e-8, overdisp);
  //     if (i == start + 1) {
  //       time_term[i] ~ normal(- autocor * sum(time_term[(start + 1):finish]), ar_scale);
  //     } else {
  //       time_term[i] ~ normal(autocor * time_term[i - 1], ar_scale);
  //     }
  //   }
  // }

}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real rate_ = exp(fmax(fmin(log_rate[n], 12.0), -12.0));
    real phi_ = fmin(fmax(overdisp, 0.01), 300.0);
    log_lik[n] = neg_binomial_2_lpmf(y[n] | rate_, phi_);
  }
}


// generated quantities {
//   vector[N] log_lik;
//   real phi_ = fmin(fmax(overdisp, 0.01), 300.0);
//   // for likelihood replace starts with - sum of all
//   for (j in 1:M) {  // for every county
//     int start = county_brks[j] + 1;
//     int finish = county_brks[j + 1];
//     // for t0 ignore the parameter and treat it as the negative sum of all others to guarantee zero-sum
//     real sum_all_but_first = sum(time_term[(start + 1):finish]);
//     real rate_ = exp(fmax(fmin(log_rate[start] - sum_all_but_first, 12.0), -12.0)) + 1e-8;
//     log_lik[start] = neg_binomial_2_lpmf(y[start] | rate_, phi_);
//     // for the rest do usual ar(1) but make sure tie it with t0 correctly
//     for (i in (start + 1):finish) {
//       real ratei_ = exp(fmax(fmin(log_rate[i] + time_term[i], 12.0), -12.0)) + 1e-8;
//       log_lik[i] = neg_binomial_2_lpmf(y[i] | ratei_, phi_);
//     }
//   }
// }
