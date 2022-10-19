data {
  int<lower=0> n_preds;
  int<lower=0> n_thr;
  int<lower=0> n_models;
  int<lower=0> n_intervals;
  vector<lower=0, upper=1>[n_thr] thresholds;
  matrix<lower=0, upper=n_preds>[n_models, n_thr] positives;
  matrix<lower=0>[n_models, n_thr] pos_post1;
  matrix<lower=0>[n_models, n_thr] pos_post2;
  // Exposure time within each interval, for each desired survival time
  row_vector[n_intervals] time_exposed;
  matrix[n_intervals, n_thr] posterior_alpha[n_models];
  matrix[n_intervals, n_thr] posterior_beta[n_models];
  vector[n_intervals] posterior_alpha0;
  vector[n_intervals] posterior_beta0;
}

parameters {
  matrix<lower=0>[n_intervals, n_thr] lambda[n_models]; 
  vector<lower=0, upper=1>[n_thr] positivity[n_models];
  vector<lower=0>[n_intervals] lambda0;
}

transformed parameters {
  vector[n_thr] St_positives[n_models];
  real St_marginal;
  
  for (model_j in 1:n_models) {
    St_positives[model_j] = exp(-(time_exposed * lambda[model_j]))'; // S(t) = exp(-H(t)), matrix * vector product
  }
  
  St_marginal = exp(-(time_exposed*lambda0))';
  
}

model {
  // one lambda parameter for each interval
  // computed using analytical posterior distribution
  // vectorized notation
  for (model_j in 1:n_models) {
    for (thr_k in 1:n_thr) {
      lambda[model_j][,thr_k] ~ gamma(posterior_alpha[model_j][,thr_k], posterior_beta[model_j][,thr_k]);
      positivity[model_j][thr_k] ~ beta(pos_post1[model_j, thr_k], pos_post2[model_j, thr_k]);
    }
  }
  
  lambda0 ~ gamma(posterior_alpha0, posterior_beta0);
  
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
    treat_all[thr_i] = (1-St_marginal) - St_marginal*odds_thr; // (1 - St_marginal)*1 - St_marginal*1*odds(thr)
    for (model_j in 1:n_models) {
      // NB for threshold i, model j
      net_benefit[thr_i, model_j] = (1-St_positives[model_j][thr_i])*positivity[model_j][thr_i] - St_positives[model_j][thr_i]*positivity[model_j][thr_i]*odds_thr;
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
