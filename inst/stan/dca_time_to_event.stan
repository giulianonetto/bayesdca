data {
  int<lower=0> n_thr;
  int<lower=0> n_models;
  int<lower=0> n_intervals;  // number of time intervals with constant hazard
  vector<lower=0, upper=1>[n_thr] thresholds;
  matrix<lower=0>[n_models, n_thr] pos_post1;  // posterior pars for P(phat > t)
  matrix<lower=0>[n_models, n_thr] pos_post2;
  row_vector[n_intervals] time_exposed;  // exposure time within each time interval
  matrix[n_intervals, n_thr] posterior_alpha[n_models]; // posterior pars for constant hazards, given phat > t
  matrix[n_intervals, n_thr] posterior_beta[n_models];
  vector[n_intervals] posterior_alpha0;  // posterior pars for constant hazards (marginal)
  vector[n_intervals] posterior_beta0;
  int<lower=1> other_models_indices[n_models, n_models-1];
}

transformed data {
  vector<lower=0>[n_thr] odds_thrs;
  for (i in 1:n_thr) {
    odds_thrs[i] = thresholds[i] / (1 - thresholds[i]);
  }
}

parameters {
  matrix<lower=0>[n_intervals, n_thr] lambda[n_models];
  vector<lower=0, upper=1>[n_thr] positivity[n_models];
  vector<lower=0>[n_intervals] lambda0;
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
  // Survival estimates
  vector<lower=0, upper=1>[n_thr] St_positives[n_models];
  real<lower=0, upper=1> St_marginal;
  // Treat all (same for all models)
  vector[n_thr] treat_all;
  // Record higher net benefit fot each threshold
  matrix[n_thr, n_models] highest_nb_other_than_model_j;
  // Net benefir calculation for each threshold, each model
  matrix[n_thr, n_models] net_benefit;
  // Delta NB calculation (against treat all/none) for each threshold, each model
  matrix[n_thr, n_models] delta;
  // Probability of being useful
  // for each threshold, each model
  matrix<lower=0, upper = 1>[n_thr, n_models] prob_better_than_soc;
  // Treat all
  St_marginal = exp(-(time_exposed*lambda0))';
  treat_all = (1-St_marginal) - St_marginal*odds_thrs; // (1 - St_marginal)*1 - St_marginal*1*odds(thr)

  for (model_j in 1:n_models) {
    // NB for model j, all thresholds
    St_positives[model_j] = exp(-(time_exposed * lambda[model_j]))';
    net_benefit[, model_j] = (1-St_positives[model_j]).*positivity[model_j] - St_positives[model_j].*positivity[model_j].*odds_thrs;
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
