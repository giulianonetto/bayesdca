#' Fit Bayesian Decision Curve Analysis using Stan for list of models or binary tests
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
#' @param prior_p,prior_se,prior_sp Prior parameters for prevalence, sensitivity, and specificity (numeric matrices of size `n_thr` by `n_models_or_tests`).
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#'
.dca_stan_list <- function(n_thr,
                           n_models_or_tests,
                           N,
                           d,
                           tp,
                           tn,
                           thresholds,
                           prior_p1, prior_p2,
                           prior_Se1, prior_Se2,
                           prior_Sp1, prior_Sp2,
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
    prior_Sp2 = prior_Sp2
  )

  .model <- stanmodels$dca_list_model
  stanfit <- rstan::sampling(.model, data = standata,
                             refresh = refresh, ...)
  return(stanfit)
}


#' Bayesian Decision Curve Analysis for a list of Predictive Models and Binary tests
#'
#' @export
#' @description Estimate decision curves for a list of predictive models and/or
#' binary tests all at once. Necessary to make comparative inferences across
#' multiple models or tests
#' using their corresponding posterior draws.
#' @param .data A data.frame with an `outcomes` column (0 or 1 for each individual)
#' and one or more columns with predicted probabilities from each of desired list
#' of predictive models, or with 0 or 1 indicator from each of desired list of
#' binary tests.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default is FALSE).
#' @param keep_draws Logical indicating whether to keep posterior
#' draws from `stanfit` object (default is TRUE).
#' @param prior_p,prior_se,prior_sp Non-negative shape values for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively.
#' Default is uniform prior for all parameters - Beta(1, 1).
#' A single vector of the form `c(a, b)` can be provided for each.
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `DCAList`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' head(PredModelData)
#' fit <- dca_list(.data = PredModelData)
#' plot(fit)
dca_list <- function(.data,
                     thresholds = seq(0.01, 0.5, 0.01),
                     keep_draws = TRUE,
                     keep_fit = FALSE,
                     prior_p = NULL,
                     prior_se = NULL,
                     prior_sp = NULL,
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
  threshold_data <- get_thr_data_list(.data = .data,
                                      thresholds = thresholds)
  priors <- get_prior_parameters(prior_p = prior_p,
                                 prior_se = prior_se,
                                 prior_sp = prior_sp,
                                 n_thresholds = n_thresholds,
                                 n_models_or_tests = n_models_or_tests)
  tp <- threshold_data %>%
    dplyr::select(thresholds, model_or_test, tp) %>%
    tidyr::pivot_wider(names_from = model_or_test,
                       values_from = tp) %>%
    dplyr::select(-thresholds) %>%
    as.matrix()

  tn <- threshold_data %>%
    dplyr::select(thresholds, model_or_test, tn) %>%
    tidyr::pivot_wider(names_from = model_or_test,
                       values_from = tn) %>%
    dplyr::select(-thresholds) %>%
    as.matrix()

  fit <- .dca_stan_list(
    n_thr = n_thresholds,
    n_models_or_tests = n_models_or_tests,
    N = rep(N, n_thresholds),
    d = rep(d, n_thresholds),
    tp = tp,
    tn = tn,
    thresholds = thresholds,
    prior_p1 = priors[['p1']],
    prior_p2 = priors[['p2']],
    prior_Se1 = priors[['Se1']],
    prior_Se2 = priors[['Se2']],
    prior_Sp1 = priors[['Sp1']],
    prior_Sp2 = priors[['Sp2']],
    refresh = refresh,
    ...
  )

  dca_summary <- extract_dca_list_summary(fit = fit,
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
    output_data[['draws']] <- extract_dca_draws(fit = fit,
                                                model_or_test_names = model_or_test_names)
  }

  .output <- structure(output_data, class = "BayesDCAList")
  return(.output)
}

#' @title Get Threshold Performance Data for DCA list
#'
#' @param outcomes Integer vector (0 or 1) with binary outcomes.
#' @param predictions Numeric vector with predicted probabilities.
#' @importFrom magrittr %>%
get_thr_data_list <- function(.data,
                              thresholds = seq(0.01, 0.5, 0.01)) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }
  outcomes <- .data[['outcomes']]
  model_or_test_names <- colnames(.data)[-1]
  N <- length(outcomes)
  d <- sum(outcomes)

  thr_data <- purrr::map(model_or_test_names, ~ {
    .predictions <- .data[[.x]]
    .thr_data <- tibble::tibble(
      model_or_test = .x,
      N = N, d = d, thresholds = thresholds
    ) %>%
      dplyr:::mutate(
        thr_perf = purrr::map(thresholds, function(.thr) {
          # for binary tests, this gives the same tp/tn for all thresholds
          tp <- sum(.predictions[outcomes == 1] >= .thr)
          tn <- sum(.predictions[outcomes == 0] < .thr)
          return(list(tp = tp, tn = tn))
        })
      ) %>%
      tidyr::unnest_wider(col = thr_perf) %>%
      dplyr::select(model_or_test, N, d, tp, tn, thresholds)
  }) %>%
    dplyr::bind_rows()

  return(thr_data)
}

#' @title Get prior parameters formated for Stan model
#'
#' @param n_thresholds Number of thresholds (int.).
#' @param n_models_or_tests Number of models or tests (int.).
#' @param prior_p,prior_se,prior_sp Non-negative shape values for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively.
#' Default is uniform prior for all parameters - Beta(1, 1).
#' A single vector of the form `c(a, b)` can be provided for each.
#' @importFrom magrittr %>%
get_prior_parameters <- function(n_thresholds,
                                 n_models_or_tests,
                                 prior_p = NULL,
                                 prior_se = NULL,
                                 prior_sp = NULL) {

  if (is.null(prior_p)) prior_p <- c(1, 1)
  if (is.null(prior_se)) prior_se <- c(1, 1)
  if (is.null(prior_sp)) prior_sp <- c(1, 1)

  stopifnot(
    length(prior_p) == 2 & is.vector(prior_p)
  )
  stopifnot(
    length(prior_se) == 2 & is.vector(prior_se)
  )
  stopifnot(
    length(prior_sp) == 2 & is.vector(prior_sp)
  )

  se1 <- sapply(1:n_models_or_tests, function(i) rep(prior_se[1], n_thresholds))
  se2 <- sapply(1:n_models_or_tests, function(i) rep(prior_se[2], n_thresholds))
  sp1 <- sapply(1:n_models_or_tests, function(i) rep(prior_sp[1], n_thresholds))
  sp2 <- sapply(1:n_models_or_tests, function(i) rep(prior_sp[2], n_thresholds))

  .priors <- list(
    p1 = prior_p[1], p2 = prior_p[2],
    Se1 = se1, Se2 = se2,
    Sp1 = sp1, Sp2 = sp2
  )
  return(.priors)
}

#' @title Get summary from BayesDCA fit
#'
#' @param fit A stanfit object.
#' @param model_or_test_names Vector of names of models or binary tests under assessment.
#' @param summary_probs Numeric vector giving probabilities for credible interval.
#' @param thresholds Vector of thresholds for DCA.
#' @importFrom magrittr %>%
extract_dca_list_summary <- function(fit,
                                     model_or_test_names,
                                     summary_probs = c(0.025, 0.975),
                                     thresholds = seq(0.01, 0.5, 0.01)) {

  # overall summary
  fit_summary <- rstan::summary(
    fit, probs = summary_probs
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::as_tibble(rownames = "par_name") %>%
    dplyr::select(-se_mean) %>%
    dplyr::rename(estimate := mean) %>%
    dplyr::mutate(
      threshold_ix = stringr::str_extract(par_name, "\\[\\d+") %>%
        stringr::str_remove("\\[") %>%
        as.integer(),
      model_or_test_ix = stringr::str_extract(par_name, ",\\d+\\]") %>%
        stringr::str_remove_all("\\]|,") %>%
        as.integer(),
      threshold = thresholds[threshold_ix],
      model_or_test_name = model_or_test_names[model_or_test_ix]
    ) %>%
    dplyr::select(par_name, threshold, model_or_test_name,
                  dplyr::everything(), -dplyr::contains("ix"))

  # prevalence
  prevalence <- fit_summary %>%
    dplyr::filter(
      par_name == "p"
    ) %>%
    dplyr::select(-threshold, -model_or_test_name)

  # sensitivity
  sensitivity <- fit_summary %>%
    dplyr::filter(
      stringr::str_detect(par_name, "Se")
    )

  # specificity
  specificity <- fit_summary %>%
    dplyr::filter(
      stringr::str_detect(par_name, "Sp")
    )

  # net benefit
  net_benefit <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "net_benefit"))

  # treat all (net benefit for treat all strategy)
  treat_all <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "treat_all")) %>%
    dplyr::select(-model_or_test_name)

  .summary <- structure(
    list(
      net_benefit = net_benefit,
      treat_all = treat_all,
      prevalence = prevalence,
      sensitivity = sensitivity,
      specificity = specificity,
      model_or_test_names = model_or_test_names
    ), class = "BayesDCASummary"
  )
  return(.summary)
}

#' @title Get posterior draws from DCA stanfit
#'
#' @param fit A stanfit object.
#' @param model_or_test_names Vector of names of models or binary tests under assessment.
extract_dca_draws <- function(fit,
                              model_or_test_names) {

  stan_draws <- rstan::extract(fit)

  .draws <- list(
    p = stan_draws$p,
    treat_all = stan_draws$treat_all,
    net_benefit = list(),
    delta_default = list(),
    prob_better_than_default = list(),
    Se = list(),
    Sp = list()
  )

  for (i in seq_along(model_or_test_names)) {
    .name <- model_or_test_names[i]
    .draws[['net_benefit']][[.name]] <- stan_draws$net_benefit[,,i]
    .draws[['delta_default']][[.name]] <- stan_draws$delta_default[,,i]
    .draws[['prob_better_than_default']][[.name]] <- stan_draws$prob_better_than_soc[,,i]
    .draws[['Se']][[.name]] <- stan_draws$Se[,,i]
    .draws[['Sp']][[.name]] <- stan_draws$Sp[,,i]
  }

  return(.draws)
}

#' @title Print BayesDCAList
#'
#' @param obj BayesDCAList object
#' @export
print.BayesDCAList <- function(obj, ...) {
  .data <- paste0(
    "N=", unique(obj$threshold_data$N), "\tD=", unique(obj$threshold_data$d)
  )
  cat(
    paste0(
      c(
        "BayesDCAList\n",
        paste0("Number of thresholds: ", length(obj$thresholds)),
        "\nRaw data:", .data, 
        "Models or tests: ",
        paste0(obj$model_or_test_names, collapse = ", ")
      ),
      collapse = "\n"
    )
  )
}
