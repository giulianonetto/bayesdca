#' Fit Bayesian Decision Curve Analysis using Stan for time-to-event models or tests
#'
#' @param n_thr Number of thresholds (int.).
#' @param n_models_or_tests Number of models or binary tests (int.).
#' @param N Sample size (vector of integers of length `n_thr`).
#' @param d Diseased: number of diseased persons or events (vector of integers of length `n_thr`).
#' @param tp True Positives: number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model (matrix of integers of size `n_thr` by `n_models_or_tests`).
#' @param tn True Negatives: number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model (matrix of integers of size `n_thr` by `n_models_or_tests`).
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0.01, 0.5, 0.01)`).
#' @param N_ext,d_ext External sample size and number of diseased individuals (or cases), respectively, used to adjust prevalence.
#' @param prior_p,prior_se,prior_sp Prior parameters for prevalence, sensitivity, and specificity (numeric matrices of size `n_thr` by `n_models_or_tests`).
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#'
.dca_stan_tte <- function(n_thr,
                          n_models_or_tests,
                          N,
                          d,
                          tp,
                          tn,
                          thresholds,
                          prior_p1, prior_p2,
                          prior_Se1, prior_Se2,
                          prior_Sp1, prior_Sp2,
                          N_ext,
                          d_ext,
                          refresh = 0, ...) {

  thresholds <- pmin(thresholds, 0.999)  # odds(1) = Inf

  standata <- list(
    n_thr = n_thr,
    n_models = n_models_or_tests,
    N = N,
    d = d,
    tp = tp,
    tn = tn,
    thresholds = thresholds,
    prior_p1 = prior_p1,
    prior_p2 = prior_p1,
    prior_Se1 = prior_Se1,
    prior_Se2 = prior_Se2,
    prior_Sp1 = prior_Sp1,
    prior_Sp2 = prior_Sp2,
    N_ext = N_ext,
    d_ext = d_ext
  )

  standata <- list(
    N_intervals = n_distinct(df_summary$j),
    N_survival_times = length(.survival_times),
    time_exposed = .time_exposed,
    posterior_alpha = .posterior$alpha,
    posterior_beta = .posterior$beta
  )

  .model <- stanmodels$dca_list_model
  stanfit <- rstan::sampling(.model, data = standata,
                             refresh = refresh, ...)
  return(stanfit)
}
