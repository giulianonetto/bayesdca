data {
  int<lower=0> N;
  int<lower=0> n_thr;                                       // number of decision thresholds
  int<lower=0> n_strategies;                                    // number of decision strategies (models or tests)
  int<lower=0> n_event_times_stacked;
  int<lower=0> n_censored_times_stacked;
  real<lower=0> event_times_stacked[n_event_times_stacked];
  real<lower=0> censored_times_stacked[n_censored_times_stacked];
  int<lower=0> event_times_start_positions[n_strategies, n_thr];
  int<lower=0> censored_times_start_positions[n_strategies, n_thr];
  int<lower=0> event_times_sizes[n_strategies, n_thr];
  int<lower=0> censored_times_sizes[n_strategies, n_thr];
  real<lower=0> prediction_time;
  vector<lower=0, upper=1>[n_thr] thresholds;               // decision thresholds
  matrix<lower=0>[n_strategies, n_thr] pos_post1;               // posterior shape1 for P(phat > t)
  matrix<lower=0>[n_strategies, n_thr] pos_post2;               // posterior shape2 for P(phat > t)
  int<lower=0> total_event_times;
  int<lower=0> total_censored_times;
  vector<lower=0>[total_event_times] event_times_marginal;
  vector<lower=0>[total_censored_times] censored_times_marginal;
  int<lower=1> other_models_indices[n_strategies, n_strategies-1];  // for each model, indices to capture all models but itself
  int<lower=1, upper=2> shape_prior;
  int<lower=1, upper=2> scale_prior;
  real<lower=0> shape_prior_pars[3];
  real<lower=0> scale_prior_pars[3];
  int<lower=0, upper=1> prior_only;
}

transformed data {
  vector<lower=0>[n_thr] odds_thrs;

  for (m in 1:n_thr) {
    odds_thrs[m] = thresholds[m] / (1 - thresholds[m]);
  }           
}

parameters {
    real<lower=0> shape[n_strategies, n_thr];
    real<lower=0> scale[n_strategies, n_thr];
    real<lower=0> shape_marginal;
    real<lower=0> scale_marginal;
}

model {
    for (j in 1:n_strategies) {
        for (i in 1:n_thr) {
            if (prior_only == 0) {
                if (event_times_sizes[j, i] > 0) {
                    segment(event_times_stacked, event_times_start_positions[j, i], event_times_sizes[j, i]) ~ weibull(shape[j, i], scale[j, i]);
                }  // else: no positives for this model at this threshold with observed event times
                if (censored_times_sizes[j, i] > 0) {
                    target += weibull_lccdf(segment(censored_times_stacked, censored_times_start_positions[j, i], censored_times_sizes[j, i]) | shape[j, i], scale[j, i] );
                }  // else: no positives for this model at this threshold with censored event times
            }
            // Prior statements
            if (shape_prior == 1) {
                shape[j, i] ~ student_t(shape_prior_pars[1], shape_prior_pars[2], shape_prior_pars[3]);
            } else {
                shape[j, i] ~ gamma(shape_prior_pars[1], shape_prior_pars[2]);
            }

            if (scale_prior == 1) {
                scale[j, i] ~ student_t(scale_prior_pars[1], scale_prior_pars[2], scale_prior_pars[3]);
            } else {
                scale[j, i] ~ gamma(scale_prior_pars[1], scale_prior_pars[2]);
            }
        }
    }

    if (prior_only == 0) {
      event_times_marginal ~ weibull(shape_marginal, scale_marginal);
      target += weibull_lccdf(censored_times_marginal | shape_marginal, scale_marginal);
    }

    if (shape_prior == 1) {
        shape_marginal ~ student_t(shape_prior_pars[1], shape_prior_pars[2], shape_prior_pars[3]);
    } else {
        shape_marginal ~ gamma(shape_prior_pars[1], shape_prior_pars[2]);
    }

    if (scale_prior == 1) {
        scale_marginal ~ student_t(scale_prior_pars[1], scale_prior_pars[2], scale_prior_pars[3]);
    } else {
        scale_marginal ~ gamma(scale_prior_pars[1], scale_prior_pars[2]);
    }

}

generated quantities {
  vector<lower=0, upper=1>[n_thr] positivity[n_strategies];
  real<lower=0, upper=1> St_marginal;
  real<lower=0, upper=1> St_positives[n_strategies, n_thr];
  // Treat all (same for all models)
  vector<upper=1>[n_thr] treat_all;
  // Record highest net benefit for each threshold among all models other than model j
  matrix[n_thr, n_strategies] best_competitor_nb;
  // Net benefir calculation for each threshold, each model
  matrix<upper=1>[n_thr, n_strategies] net_benefit;
  // Delta NB calculation (against treat all/none) for each threshold, each model
  matrix[n_thr, n_strategies] delta_best;
  matrix[n_thr, n_strategies] delta_useful;
  // Probability of being best/useful for each threshold, each model
  matrix<lower=0, upper=1>[n_thr, n_strategies] p_best;
  matrix<lower=0, upper=1>[n_thr, n_strategies] p_useful;
  // Treat all
  St_marginal = exp((-1) * pow( prediction_time / scale_marginal, shape_marginal));
  treat_all = (1-St_marginal) * 1 - St_marginal * 1 * odds_thrs;

  for (j in 1:n_strategies) {
    for (i in 1:n_thr) {
        positivity[j][i] = beta_rng(pos_post1[j, i], pos_post2[j, i]);    
        St_positives[j, i] = exp( (-1) * pow(prediction_time / scale[j, i], shape[j, i]) );  // survival at prediction time
        net_benefit[i, j] = ((1 - St_positives[j, i]) * positivity[j, i]) - ((St_positives[j, i] * positivity[j, i]) * odds_thrs[i]);
    }
  }

  for (i in 1:n_thr) {
    real best_among_treat_all_or_none = fmax(0, treat_all[i]);

    for (j in 1:n_strategies) {

      if (n_strategies > 1) {
        real best_among_other_models = max(net_benefit[i, other_models_indices[j]]);
        best_competitor_nb[i, j] = fmax(best_among_treat_all_or_none, best_among_other_models);

      } else {
        best_competitor_nb[i, j] = best_among_treat_all_or_none;
      }

      // P(best)
      p_best[i, j] = net_benefit[i, j] > best_competitor_nb[i, j];
      // P(useful)
      p_useful[i, j] = net_benefit[i, j] > best_among_treat_all_or_none;
      // Delta against best strategy other than j^th strategy
      delta_best[i, j] = net_benefit[i, j] - best_competitor_nb[i, j];
      // Delta agains Treat all/none
      delta_useful[i, j] = net_benefit[i, j] - best_among_treat_all_or_none;
    }
  }
}
