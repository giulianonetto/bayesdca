#' Fit Bayesian Decision Curve Analysis using Stan for survival outcomes
#'
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#'
.dca_stan_surv <- function(n_thr,
                           n_models_or_tests,
                           n_intervals,
                           thresholds,
                           time_exposed,
                           posterior_alpha,
                           posterior_beta,
                           posterior_alpha0,
                           posterior_beta0,
                           pos_post1,
                           pos_post2,
                           refresh = 0, ...) {

  thresholds <- pmin(thresholds, 0.9999)  # odds(1) = Inf

  standata <- list(
    n_thr = n_thr,
    n_models = n_models_or_tests,
    n_intervals = n_intervals,
    thresholds = thresholds,
    time_exposed = time_exposed,
    posterior_alpha = posterior_alpha,
    posterior_beta = posterior_beta,
    posterior_alpha0 = posterior_alpha0,
    posterior_beta0 = posterior_beta0,
    pos_post1 = pos_post1,
    pos_post2 = pos_post2,
    refresh = refresh
  )

  .model <- stanmodels$dca_time_to_event
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
                     prediction_time,
                     thresholds = seq(0.01, 0.5, 0.01),
                     keep_draws = TRUE,
                     keep_fit = FALSE,
                     summary_probs = c(0.025, 0.975),
                     cutpoints = NULL,
                     refresh = 0,
                     ...) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }

  stopifnot(
    "'outcomes' column must be Surv object. " = survival::is.Surv(.data[['outcomes']])
  )

  # avoid thresholds in {0, 1}
  thresholds <- thresholds %>%
    pmin(0.99999) %>%
    pmax(0.00001)

  # preprocess .data
  prediction_data <- data.frame(.data[, -1])
  surv_data <- data.frame(
    .time = unname(.data[['outcomes']][,1]),   # observed time-to-event
    .status = unname(.data[['outcomes']][,2])  # 1 if event, 0 if censored
  )

  event_times <- surv_data$.time[surv_data$.status > 0]
  model_or_test_names <- colnames(prediction_data)

  # preprocess cutpoints
  if (is.null(cutpoints)) {
    cutpoints <- quantile(
      event_times, probs = c(0.1, 0.5, 0.9)
    )
  }

  # make sure zero is included and cutpoints are correctly ordered
  cutpoints <- sort(unique(c(0, cutpoints)))

  stopifnot(
    "prediction_time must be less than max(cutpoints)" = prediction_time < max(cutpoints)
  )


  events_per_interval <- table(cut(event_times, c(cutpoints, Inf), include.lowest = T))
  if (!(all(events_per_interval >= 5))) {
    msg <- paste0(
      "Please updade cutpoints, too few events per interval (at least five needed):\n",
      paste(names(events_per_interval), events_per_interval, sep = ": ", collapse = "  ")
    )
    stop(msg)
  }

  n_models_or_tests <- ncol(prediction_data)
  n_thresholds <- length(thresholds)
  n_intervals = length(cutpoints)
  time_exposed <- get_survival_time_exposed(prediction_time, cutpoints)

  posterior_surv_pars <- get_survival_posterior_parameters(
    .prediction_data = prediction_data,
    .surv_data = surv_data,
    .models_or_tests = colnames(prediction_data),
    .cutpoints = cutpoints,
    .thresholds = thresholds
  )
  posterior_surv_pars0 <- get_survival_posterior_parameters(
    .prediction_data = prediction_data,
    .surv_data = surv_data,
    .models_or_tests = colnames(prediction_data),
    .cutpoints = cutpoints,
    .thresholds = 0
  )
  posterior_positivity_pars <- get_positivity_posterior_parameters(
    .prediction_data = prediction_data,
    .thresholds = thresholds
  )

  fit <- .dca_stan_surv(
    n_thr = n_thresholds,
    n_models = n_models_or_tests,
    n_intervals = n_intervals,
    thresholds = thresholds,
    time_exposed = time_exposed, # get_time_exposed(pred_time, .cutpoints) %>% as.vector(),
    posterior_alpha = posterior_surv_pars$.alpha, #list(sapply(1:length(thr), \(i) posterior_pars[[i]]$alpha)),
    posterior_beta = posterior_surv_pars$.beta, #list(sapply(1:length(thr), \(i) posterior_pars[[i]]$beta)),
    posterior_alpha0 = posterior_surv_pars0$.alpha, #posterior_pars0[[1]]$alpha,
    posterior_beta0 = posterior_surv_pars0$.beta,#posterior_pars0[[1]]$beta,
    pos_post1 = posterior_positivity_pars$.shape1, #0.5 + t(sapply(thr, \(i) sum(preds>i))),
    pos_post2 = posterior_positivity_pars$.shape2, #0.5 + length(preds) - t(sapply(thr, \(i) sum(preds>i))),
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
    event_times = event_times,
    censor_times = surv_data$.time[surv_data$.status == 0],
    .data = .data,
    cutpoints = cutpoints,
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
