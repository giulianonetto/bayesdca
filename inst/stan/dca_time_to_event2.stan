data {
  int<lower=0> n_thr;                                       // number of decision thresholds
  int<lower=0> n_models;                                    // number of decision strategies (models or tests)
  int<lower=0> n_intervals;                                 // number of time intervals with constant hazard
  vector<lower=0, upper=1>[n_thr] thresholds;               // decision thresholds
  matrix<lower=0>[n_intervals, n_thr] time_exposed[n_models];         // t_ik's for each interval
  row_vector[n_intervals] time_exposed_prediction;  // exposure time within each time interval for prediciton of survival
  int<lower=0> deaths[n_models, n_intervals, n_thr];                 // d_ik's for each model
  matrix<lower=0>[n_models, n_thr] pos_post1;               // posterior shape1 for P(phat > t)
  matrix<lower=0>[n_models, n_thr] pos_post2;               // posterior shape2 for P(phat > t)
  vector<lower=0>[n_intervals] posterior_alpha0;                     // posterior alpha for constant hazards (marginal, unsed in Treat all)
  vector<lower=0>[n_intervals] posterior_beta0;                      // posterior beta for constant hazards (marginal, unsed in Treat all)
  int<lower=1> other_models_indices[n_models, n_models-1];  // for each model, indices to capture all models but itself
  int<lower=0, upper=1> prior_only;
}

transformed data {
  vector<lower=0>[n_thr] odds_thrs;
  for (i in 1:n_thr) {
    odds_thrs[i] = thresholds[i] / (1 - thresholds[i]);
  }
}

parameters {
  matrix<lower=0>[n_intervals, n_thr] lambda[n_models];     // interval- and threshold- specific hazard for each model
  vector<lower=0>[n_thr] sd_hazard[n_models];                      // par2 for hierarchical prior on lambda
  vector<lower=0, upper=1>[n_thr] positivity[n_models];     // prob. of positive prediction at each threshold for each model
  vector<lower=0>[n_intervals] lambda0;                     // used for treat all strategy
}

transformed parameters {
  vector<lower=0, upper=1>[n_thr] St_positives[n_models];   // conditional survival given positive prediction (prediction above decision threshold)
  real<lower=0, upper=1> St_marginal;

  for (model_j in 1:n_models) {
    St_positives[model_j] = exp(-(time_exposed_prediction * lambda[model_j]))'; // S(t) = exp(-H(t)), matrix * vector product
  }

  St_marginal = exp(-(time_exposed_prediction * lambda0))';

}

model {
  for (model_j in 1:n_models) {
    sd_hazard[model_j] ~ normal(0, 0.15);
    for (thr_m in 1:n_thr) {
      // P(+ | threshold) varies by model and threshold only
      positivity[model_j][thr_m] ~ beta(pos_post1[model_j, thr_m], pos_post2[model_j, thr_m]);
      for (k in 1:n_intervals) {
        lambda[model_j][k, thr_m] ~ exponential(1); // prior on hazards
        if (prior_only == 0) {
            deaths[model_j, k, thr_m] ~ poisson(time_exposed[model_j][k, thr_m] * lambda[model_j][k, thr_m]); // likelihood
        }
      }
    }
  }
  lambda0 ~ gamma(posterior_alpha0, posterior_beta0);
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
  // Probability of being useful
  // for each threshold, each model
  matrix<lower=0, upper = 1>[n_thr, n_models] prob_better_than_soc;
  // Treat all
  treat_all = (1-St_marginal) - St_marginal*odds_thrs; // (1 - St_marginal)*1 - St_marginal*1*odds(thr)

  for (model_j in 1:n_models) {
    // NB for model j, all thresholds
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
