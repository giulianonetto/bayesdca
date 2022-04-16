data {
  // Input Data
  int<lower = 1> n_thr;  // number of thresholds
  int<lower = 0> n_models;  // number of models (or binary tests) to evaluate
  vector<lower = 0>[n_thr] N;  // sample size (same for each threshold)
  matrix<lower = 0>[n_thr, n_models] tp; // true positives for each threshold, each model
  matrix<lower = 0>[n_thr, n_models] tn; // true negatives for each threshold, each model
  vector<lower = 0>[n_thr] d;  // diseased persons or outcome cases (same for each threshold)
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
  // Analytical posterior parameters for prevalence
  post_p1 = d[1] + prior_p1;
  post_p2 = N[1] - d[1] + prior_p2;
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
  // Net benefir calculation for each threshold, each model
  matrix[n_thr, n_models] net_benefit;
  // Delta NB calculation (against treat all/none) for each threshold, each model
  matrix[n_thr, n_models] delta;
  // Probability that delta NB is greater than zero - model better than treat all/none ("Standard of Care")
  // for each threshold, each model
  matrix<lower=0, upper = 1>[n_thr, n_models] prob_better_than_soc;

  for (thr_i in 1:n_thr) {
    // tmp variable for odds(threshold)
    real odds_thr = thresholds[thr_i]/(1-thresholds[thr_i]);
    // Treat all
    treat_all[thr_i] = 1*p - (1-p)*(1-0)*odds_thr;
    for (model_j in 1:n_models) {
      // NB for threshold i, model j
      net_benefit[thr_i, model_j] = Se[thr_i, model_j]*p - (1-p)*(1-Sp[thr_i, model_j])*odds_thr;
      // Delta NB for threshold i, model j:
      // against treat all if treat all is positive
      // against treat none otherwise
      if (treat_all[thr_i] > 0) {
        delta[thr_i, model_j] = net_benefit[thr_i, model_j] - treat_all[thr_i];
      } else {
        delta[thr_i, model_j] = net_benefit[thr_i, model_j];  // "minus treat none" = "minus zero"
      }
      // P(delta NB > 0) = probability model is better than treat all/none
      prob_better_than_soc[thr_i, model_j] = delta[thr_i, model_j] > 0;
    }
  }
}
