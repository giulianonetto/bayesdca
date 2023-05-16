data {
  // Input Data
  int<lower = 1> n_thr;  // number of thresholds
  int<lower = 0> n_models;  // number of models (or binary tests) to evaluate
  vector<lower = 0>[n_thr] N;  // sample size (same for each threshold)
  matrix<lower = 0>[n_thr, n_models] tp; // true positives for each threshold, each model
  matrix<lower = 0>[n_thr, n_models] tn; // true negatives for each threshold, each model
  vector<lower = 0>[n_thr] d;  // diseased persons or outcome cases (same for each threshold)
  vector<lower = 0>[n_thr] d_ext;  // diseased persons or outcome cases from external dataset used for prevalence adjustment (same for each threshold)
  vector<lower = 0>[n_thr] N_ext;  // sample size from external dataset used for prevalence adjustment (same for each threshold)
  vector<lower = 0, upper = 1>[n_thr] thresholds;
  // Prior parameters for prevalence
  real<lower = 0> prior_p1;
  real<lower = 0> prior_p2;
  // Prior parameters for Sensitivity (same for each threshold, all models)
  matrix<lower = 0>[n_thr, n_models] prior_Se1;
  matrix<lower = 0>[n_thr, n_models] prior_Se2;
  // Prior parameters for Specificity (same for each threshold, all models)
  matrix<lower = 0>[n_thr, n_models] prior_Sp1;
  matrix<lower = 0>[n_thr, n_models] prior_Sp2;
  int<lower=1> other_models_indices[n_models, n_models-1];
}
transformed data {
  // Posterior parameters for prevalence
  real<lower = 0> post_p1;
  real<lower = 0> post_p2;
  // Posterior parameters for Sensitivity for each threshold, each model
  matrix<lower = 0>[n_thr, n_models] post_Se1;
  matrix<lower = 0>[n_thr, n_models] post_Se2;
  // Posterior parameters for Specificity for each threshold, each model
  matrix<lower = 0>[n_thr, n_models] post_Sp1;
  matrix<lower = 0>[n_thr, n_models] post_Sp2;
  vector<lower=0>[n_thr] odds_thrs;
  // just a vector with odds of the thresholds
  for (i in 1:n_thr) {
    odds_thrs[i] = thresholds[i] / (1 - thresholds[i]);
  }
  // Analytical posterior parameters for prevalence
  if (N_ext[1] >= 1) { // it means there's external information for prevalence
    post_p1 = d_ext[1] + prior_p1;
    post_p2 = N_ext[1] - d_ext[1] + prior_p2;
  } else {
    post_p1 = d[1] + prior_p1;
    post_p2 = N[1] - d[1] + prior_p2;
  }
  // Analytical posterior parameters for Sens/Spec at threshold i, model j
  for (thr_i in 1:n_thr) {
    for (model_j in 1:n_models) {
      // Sensitivity parameters
      post_Se1[thr_i, model_j] = tp[thr_i, model_j] + prior_Se1[thr_i, model_j];
      post_Se2[thr_i, model_j] = d[1] - tp[thr_i, model_j] + prior_Se2[thr_i, model_j];
      // Specificity parameters
      post_Sp1[thr_i, model_j] = tn[thr_i, model_j] + prior_Sp1[thr_i, model_j];
      post_Sp2[thr_i, model_j] = N[1] - d[1] - tn[thr_i, model_j] + prior_Sp2[thr_i, model_j];
    }
  }
}
parameters {
  real<lower=0, upper = 1> p;
  matrix<lower=0, upper = 1>[n_thr, n_models] Se;
  matrix<lower=0, upper = 1>[n_thr, n_models] Sp;
}

transformed parameters {
  vector<lower = 0, upper = 1>[n_thr] p_vec;
  for (i in 1:n_thr) {
    p_vec[i] = p;
  }
}
model {
  // Posterior for prevalence (or outcome proportion) given parameters computed analytically
  p ~ beta(post_p1, post_p2);
  // Posterior for Sens/Spec given parameters computed analytically
  for (thr_i in 1:n_thr) {
    for (model_j in 1:n_models) {
      // Sensitivity
      Se[thr_i, model_j] ~ beta(post_Se1[thr_i, model_j], post_Se2[thr_i, model_j]);
      // Specificity
      Sp[thr_i, model_j] ~ beta(post_Sp1[thr_i, model_j], post_Sp2[thr_i, model_j]);
    }
  }
}
generated quantities {

  // Treat all (same for all models)
  vector[n_thr] treat_all;
  // Record higher net benefit fot each threshold
  matrix[n_thr, n_models] highest_nb_other_than_model_j;
  // Net benefir calculation for each threshold, each model
  matrix[n_thr, n_models] net_benefit;
  // Delta NB calculation (against treat all/none) for each threshold, each model
  matrix[n_thr, n_models] delta;
  // Probability that delta NB is greater than zero - model better than treat all/none ("Standard of Care")
  // for each threshold, each model
  matrix<lower=0, upper = 1>[n_thr, n_models] prob_better_than_soc;

  // Treat all
  treat_all = 1*p - (1-p)*(1-0)*odds_thrs;

  for (model_j in 1:n_models) {
    // NB for threshold i, model j
      net_benefit[, model_j] = Se[, model_j]*p - (1-p)*(1-Sp[, model_j]).*odds_thrs;
  }

  for (thr_i in 1:n_thr) {
    real best_among_treat_all_or_none = fmax(0, treat_all[thr_i]);

    for (model_j in 1:n_models) {

      if (n_models > 1) {
        real best_among_other_models = max(net_benefit[thr_i, other_models_indices[model_j]]);
        highest_nb_other_than_model_j[thr_i, model_j] = fmax(best_among_treat_all_or_none, best_among_other_models);

      } else {
        highest_nb_other_than_model_j[thr_i, model_j] = best_among_treat_all_or_none;
      }

      // P(useful)
      prob_better_than_soc[thr_i, model_j] = net_benefit[thr_i, model_j] > highest_nb_other_than_model_j[thr_i, model_j];
      // Delta against best strategy
      delta[thr_i, model_j] = net_benefit[thr_i, model_j] - highest_nb_other_than_model_j[thr_i, model_j];
    }
  }
}
