functions {
  // Defines the log hazard - borrowed and adapted from the survHEhmc R package :)
  real log_h(real t, real shape, real scale) {
    real log_hazard;
    log_hazard = log(shape)+(shape-1)*log(t / scale)-log(scale);
    return log_hazard;
  }
  
  // Defines the log survival
  real log_S(real t, real shape, real scale) {
    real log_survival;
    log_survival = -pow((t/scale), shape);
    return log_survival;
  }

  real censored_weibull_likelihood(real t, real d, real shape, real scale) {
    real log_lik = d * log_h(t, shape, scale) + log_S(t, shape, scale);
    return log_lik;
  }
  
  // Defines the sampling distribution
  real likelihood_lpdf(data matrix dataset, real[] shape, real[] scale, real positivity) {
    real log_lik = 0;
    int n = rows(dataset);
    int total_positives = 0;
    for (i in 1:n) {
        real time = dataset[i, 1];
        real death = dataset[i, 2];
        int is_positive = dataset[i, 3] > 0 ? 1 : 0;
        int ix = is_positive == 0 ? 1 : 2;
        log_lik += death * log_h(time, shape[ix], scale[ix]) + log_S(time, shape[ix], scale[ix]);
        if (positivity < 1) {
            log_lik += bernoulli_lpmf(is_positive | positivity);
        }        
    }
    return log_lik;
  }

}

data {
  int<lower=0> N;
  int<lower=0> n_thr;                                       // number of decision thresholds
  int<lower=0> n_models;
  vector<lower=0>[N] times;
  vector<lower=0, upper=1>[N] events;
  vector<lower=0, upper=1>[N] predictions[n_models, n_thr];    // 1 if predicted risk is above threshold, 0 otherwise
  real<lower=0> prediction_time;
  vector<lower=0, upper=1>[n_thr] thresholds;
  int<lower=1> other_models_indices[n_models, n_models-1];  // for each model, indices to capture all models but itself
  int<lower=1> shape_prior;
  int<lower=1> scale_prior;
  real<lower=0> shape_prior_pars[3];
  real<lower=0> scale_prior_pars[3];
  real<lower=0> positivity_prior_pars[2];
  int<lower=0, upper=1> prior_only;
}

transformed data {
  vector<lower=0>[n_thr] odds_thrs;
  matrix[N, 3] nb_data[n_models, n_thr];
  matrix[N, 3] nb_data_marginal;

  for (thr_m in 1:n_thr) {
    odds_thrs[thr_m] = thresholds[thr_m] / (1 - thresholds[thr_m]);
    for (model_j in 1:n_models) {
        nb_data[model_j, thr_m, , 1] = times;
        nb_data[model_j, thr_m, , 2] = events;
        nb_data[model_j, thr_m, , 3] = predictions[model_j, thr_m];
    }
  }

  nb_data_marginal[, 1] = times;
  nb_data_marginal[, 2] = events;
  nb_data_marginal[, 3] = rep_vector(1, N);
}

parameters {
    real<lower=0> shape[n_models, n_thr, 2];
    real<lower=0> scale[n_models, n_thr, 2];
    real<lower=0> shape_marginal[2];
    real<lower=0> scale_marginal[2];
    vector<lower=0, upper=1>[n_thr] positivity[n_models];
}

model {
    for (model_j in 1:n_models) {
        for (thr_m in 1:n_thr) {
            if (prior_only == 0) {
                
                target += likelihood_lpdf(nb_data[model_j, thr_m] | shape[model_j, thr_m], scale[model_j, thr_m], positivity[model_j, thr_m]);
            }
            // Prior statements
            if (shape_prior == 1) {
                target += gamma_lpdf(to_vector(shape[model_j, thr_m]) | shape_prior_pars[1], shape_prior_pars[2]);
            } else {
                target += student_t_lpdf(to_vector(shape[model_j, thr_m]) | shape_prior_pars[1], shape_prior_pars[2], shape_prior_pars[3]);
            }

            if (scale_prior == 1) {
                target += gamma_lpdf(to_vector(scale[model_j, thr_m]) | scale_prior_pars[1], scale_prior_pars[2]);
            } else {
                target += student_t_lpdf(to_vector(scale[model_j, thr_m]) | scale_prior_pars[1], scale_prior_pars[2], scale_prior_pars[3]);
            }

            target += beta_lpdf(positivity[model_j][thr_m] | positivity_prior_pars[1], positivity_prior_pars[2]);
        }
    }

    // marginal model (for Treat all, so no positivity parameter as it is exactly 1).
    if (prior_only == 0) {
      target += likelihood_lpdf(nb_data_marginal | shape_marginal, scale_marginal, 1.0);
    }

    if (shape_prior == 1) {
        target += gamma_lpdf(to_vector(shape_marginal) | shape_prior_pars[1], shape_prior_pars[2]);
    } else {
        target += student_t_lpdf(to_vector(shape_marginal) | shape_prior_pars[1], shape_prior_pars[2], shape_prior_pars[3]);
    }

    if (scale_prior == 1) {
        target += gamma_lpdf(to_vector(scale_marginal) | scale_prior_pars[1], scale_prior_pars[2]);
    } else {
        target += student_t_lpdf(to_vector(scale_marginal) | scale_prior_pars[1], scale_prior_pars[2], scale_prior_pars[3]);
    }
}

generated quantities {
  real<lower=0, upper=1> St_marginal;
  real<lower=0, upper=1> St_positives[n_models, n_thr];
  // Treat all (same for all models)
  vector<upper=1>[n_thr] treat_all;
  // Record higher net benefit fot each threshold
  matrix[n_thr, n_models] highest_nb_other_than_model_j;
  // Net benefir calculation for each threshold, each model
  matrix<upper=1>[n_thr, n_models] net_benefit;
  // Delta NB calculation (against treat all/none) for each threshold, each model
  matrix[n_thr, n_models] delta;
  // Probability of being useful
  // for each threshold, each model
  matrix<lower=0, upper=1>[n_thr, n_models] prob_better_than_soc;
  // Treat all
  St_marginal = exp(log_S(prediction_time, shape_marginal[2], scale_marginal[2]));
  treat_all = (1-St_marginal) * 1 - St_marginal * 1 * odds_thrs;

  for (model_j in 1:n_models) {
    for (thr_m in 1:n_thr) {
        real s = exp(log_S(prediction_time, shape[model_j, thr_m, 2], scale[model_j, thr_m, 2]));
        real p = positivity[model_j, thr_m];
        St_positives[model_j, thr_m] = s;
        net_benefit[thr_m, model_j] = ( (1 - s) * p ) - ( (s * p) * odds_thrs[thr_m] );
    }
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

      // P(useful or best)
      prob_better_than_soc[thr_i, model_j] = net_benefit[thr_i, model_j] > highest_nb_other_than_model_j[thr_i, model_j];
      // Delta against best strategy
      delta[thr_i, model_j] = net_benefit[thr_i, model_j] - highest_nb_other_than_model_j[thr_i, model_j];
    }
  }
}
