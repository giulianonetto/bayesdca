#' Fit Bayesian Decision Curve Analysis using Stan for survival outcomes
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
.dca_stan_surv <- function(n_preds,
                           n_thr,
                           n_models_or_tests,
                           n_intervals,
                           thresholds,
                           positives,
                           time_exposed,
                           posterior_alpha,
                           posterior_beta,
                           posterior_alpha0,
                           posterior_beta0,
                           pos_post1,
                           pos_post2,
                           refresh = 0, ...) {
  
  thresholds <- pmin(thresholds, 0.999)  # odds(1) = Inf
  
  standata <- list(
    n_preds = length(preds),
    n_thr = length(thr),
    n_models = n_models_or_tests,
    n_intervals = length(.cutpoints),
    thresholds = thresholds,
    positives = t(sapply(thr, \(i) sum(preds>i))),
    time_exposed = get_time_exposed(pred_time, .cutpoints) %>% as.vector(),
    posterior_alpha = list(sapply(1:length(thr), \(i) posterior_pars[[i]]$alpha)),
    posterior_beta = list(sapply(1:length(thr), \(i) posterior_pars[[i]]$beta)),
    posterior_alpha0 = posterior_pars0[[1]]$alpha,
    posterior_beta0 = posterior_pars0[[1]]$beta,
    pos_post1 = 0.5 + t(sapply(thr, \(i) sum(preds>i))),
    pos_post2 = 0.5 + length(preds) - t(sapply(thr, \(i) sum(preds>i)))
  )
  
  .model <- stanmodels$dca_list_model
  stanfit <- rstan::sampling(.model, data = standata,
                             refresh = refresh, ...)
  return(stanfit)
}


#' Bayesian Decision Curve Analysis for Predictive Models and Binary tests
#'
#' @export
#' @return An object of class `BayesDCASurv`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot(fit)
dca_surv <- function(.data,
                     thresholds = seq(0.01, 0.5, 0.01),
                     keep_draws = TRUE,
                     keep_fit = FALSE,
                     summary_probs = c(0.025, 0.975),
                     refresh = 0,
                     ...) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }
  model_or_test_names <- colnames(.data)[-1]
  N <- nrow(.data)
  d <- sum(.data[['outcomes']])
  n_thresholds = length(thresholds)
  n_models_or_tests = length(model_or_test_names)
  threshold_data <- .get_thr_data_list(.data = .data,
                                       thresholds = thresholds)
  priors <- .get_prior_parameters(prior_p = prior_p,
                                  prior_se = prior_se,
                                  prior_sp = prior_sp,
                                  n_thresholds = n_thresholds,
                                  n_models_or_tests = n_models_or_tests)
  time_exposed <- get_time_exposed(pred_time, .cutpoints)
  posterior_pars0 <- get_survival_baseline_posterior_parameters()
  posterior_pars <- get_survival_posterior_parameters()
  positivity_pars <- get_positivity_posterior_parameters()
  
  fit <- .dca_stan_surv(
    n_thr = n_thresholds,
    n_models = n_models_or_tests,
    n_intervals = length(.cutpoints),
    thresholds = thresholds,
    time_exposed = time_exposed, # get_time_exposed(pred_time, .cutpoints) %>% as.vector(),
    posterior_alpha = posterior_pars$alpha, #list(sapply(1:length(thr), \(i) posterior_pars[[i]]$alpha)),
    posterior_beta = posterior_pars$beta, #list(sapply(1:length(thr), \(i) posterior_pars[[i]]$beta)),
    posterior_alpha0 = posterior_pars0$alpha, #posterior_pars0[[1]]$alpha,
    posterior_beta0 = posterior_pars0$beta,#posterior_pars0[[1]]$beta,
    pos_post1 = positivity_pars$shape1, #0.5 + t(sapply(thr, \(i) sum(preds>i))),
    pos_post2 = positivity_pars$shape1, #0.5 + length(preds) - t(sapply(thr, \(i) sum(preds>i))),
    refresh = refresh,
    ...
  )
  
  dca_summary <- .extract_dca_surv_summary(fit = fit,
                                           summary_probs = summary_probs,
                                           thresholds = thresholds,
                                           model_or_test_names = model_or_test_names)
  
  output_data <- list(
    summary = dca_summary,
    thresholds = thresholds,
    .data = .data,
    threshold_data = threshold_data,
    priors = priors,
    model_or_test_names = model_or_test_names
  )
  
  if (isTRUE(keep_fit)) {
    output_data[['fit']] <- fit
  }
  
  if (isTRUE(keep_draws)) {
    output_data[['draws']] <- .extract_dca_surv_draws(fit = fit,
                                                      model_or_test_names = model_or_test_names)
  }
  
  .output <- structure(output_data, class = "BayesDCASurv")
  return(.output)
}
