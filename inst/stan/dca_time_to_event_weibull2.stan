data {
  int<lower=0> N;
  int<lower=0> n_thr;                                       // number of decision thresholds
  int<lower=0> n_models;                                    // number of decision strategies (models or tests)
  int<lower=0> n_event_times_stacked;
  int<lower=0> n_censored_times_stacked;
  real<lower=0> event_times_stacked[n_event_times_stacked];
  real<lower=0> censored_times_stacked[n_censored_times_stacked];
  int<lower=0> event_times_start_positions[n_models, n_thr];
  int<lower=0> censored_times_start_positions[n_models, n_thr];
  int<lower=0> event_times_sizes[n_models, n_thr];
  int<lower=0> censored_times_sizes[n_models, n_thr];
  real<lower=0> prediction_time;
  vector<lower=0, upper=1>[n_thr] thresholds;               // decision thresholds
  matrix<lower=0>[n_models, n_thr] pos_post1;               // posterior shape1 for P(phat > t)
  matrix<lower=0>[n_models, n_thr] pos_post2;               // posterior shape2 for P(phat > t)
  int<lower=0> total_event_times;
  int<lower=0> total_censored_times;
  vector<lower=0>[total_event_times] event_times_marginal;
  vector<lower=0>[total_censored_times] censored_times_marginal;
  int<lower=1> other_models_indices[n_models, n_models-1];  // for each model, indices to capture all models but itself
  real<lower=0> prior_sd_alpha;
  real<lower=0> prior_sd_sigma;
  int<lower=0, upper=1> prior_only;
}

transformed data {
  vector<lower=0>[n_thr] odds_thrs;

  for (m in 1:n_thr) {
    odds_thrs[m] = thresholds[m] / (1 - thresholds[m]);
  }           
}

parameters {
    real<lower=0> alpha[n_models, n_thr];
    real<lower=0> sigma[n_models, n_thr];
    real<lower=0> alpha_marginal;
    real<lower=0> sigma_marginal;
    vector<lower=0, upper=1>[n_thr] positivity[n_models];
}

model {
    for (model_j in 1:n_models) {
        for (thr_m in 1:n_thr) {
            if (prior_only == 0) {
                segment(event_times_stacked, event_times_start_positions[model_j, thr_m], event_times_sizes[model_j, thr_m]) ~ weibull(alpha[model_j, thr_m], sigma[model_j, thr_m]);
                target += weibull_lccdf(segment(censored_times_stacked, censored_times_start_positions[model_j, thr_m], censored_times_sizes[model_j, thr_m]) | alpha[model_j, thr_m], sigma[model_j, thr_m] );
            }
            // Prior statements
            alpha[model_j, thr_m] ~ student_t(3, 0, prior_sd_alpha);
            sigma[model_j, thr_m] ~ student_t(3, 0, prior_sd_sigma);
            // P(+ | prediction > threshold)
            positivity[model_j][thr_m] ~ beta(pos_post1[model_j, thr_m], pos_post2[model_j, thr_m]);    
        }
    }

    if (prior_only == 0) {
      event_times_marginal ~ weibull(alpha_marginal, sigma_marginal );
      target += weibull_lccdf(censored_times_marginal | alpha_marginal, sigma_marginal );
    }

    alpha_marginal ~ student_t(3, 0, prior_sd_alpha);
    sigma_marginal ~ student_t(3, 0, prior_sd_sigma);
    


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
  St_marginal = exp((-1) * ( prediction_time / sigma_marginal )^alpha_marginal);
  treat_all = (1-St_marginal) * 1 - St_marginal * 1 * odds_thrs;

  for (model_j in 1:n_models) {
    for (thr_m in 1:n_thr) {
        St_positives[model_j, thr_m] = exp( (-1) * ( ( prediction_time / sigma[model_j, thr_m] )^alpha[model_j, thr_m] ) );  // survival at prediction time
        net_benefit[thr_m, model_j] = ((1 - St_positives[model_j, thr_m]) * positivity[model_j, thr_m]) - ((St_positives[model_j, thr_m] * positivity[model_j, thr_m]) * odds_thrs[thr_m]);
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

      // P(useful)
      prob_better_than_soc[thr_i, model_j] = net_benefit[thr_i, model_j] > highest_nb_other_than_model_j[thr_i, model_j];
      // Delta against best strategy
      delta[thr_i, model_j] = net_benefit[thr_i, model_j] - highest_nb_other_than_model_j[thr_i, model_j];
    }
  }
}
