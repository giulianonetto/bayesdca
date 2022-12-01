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
#' @param N_ext,d_ext External sample size and number of diseased individuals (or cases), respectively, used to adjust prevalence.
#' @param prior_p,prior_se,prior_sp Prior parameters for prevalence, sensitivity, and specificity (numeric matrices of size `n_thr` by `n_models_or_tests`).
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#' @keywords internal
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
                           N_ext,
                           d_ext,
                           other_models_indices,
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
    d_ext = d_ext,
    other_models_indices = other_models_indices
  )

  .model <- stanmodels$dca_list_model
  stanfit <- rstan::sampling(.model, data = standata,
                             refresh = refresh, ...)
  return(stanfit)
}


#' Bayesian Decision Curve Analysis for Predictive Models and Binary tests
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
#' the net benefit should be computed (default is `seq(0.01, 0.5, 0.02)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default is FALSE).
#' @param keep_draws Logical indicating whether to keep posterior
#' draws from `stanfit` object (default is TRUE).
#' @param prior_p,prior_se,prior_sp Non-negative shape values for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively.
#' Default is uniform prior for all parameters - Beta(1, 1).
#' A single vector of the form `c(a, b)` can be provided for each.
#' @param priors A list with threshold- and model-specific priors
#' should contain a vector for shape1 of prevalence (named `p1`)
#' and shape2 (named `p2`). Similarly for Se1/Se2 and Sp1/Sp2, except
#' these should be matrices with as many rows as thresholds and
#' as many columns as models or tests.
#' @param constant_prior If TRUE (default), it will set a single
#' prior for all models or tests in all thresholds.
#' If FALSE, the prior will be threshold and, potentially, model/test-specific.
#' @param prior_only If set to TRUE, will produce prior DCA.
#' @param min_prior_mean,max_prior_mean Minimum
#' and maximum prior mean for sensitivity and specificity.
#' Only used if `constant_prior = FALSE`.
#' @param summary_probs Probabilities used to compute credible intervals (defaults to a 95% Cr.I.).
#' @param external_prevalence_data Vector with two positive integers giving number of diseased and
#' non-diseased individuals, respectively, from external data (e.g., if analyzing nested case-control data,
#' this is the number of cases and non-cases in the source population).
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `BayesDCAList`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot(fit)
dca <- function(.data,
                thresholds = seq(0.01, 0.5, 0.01),
                keep_draws = TRUE,
                keep_fit = FALSE,
                prior_p = NULL,
                prior_se = NULL,
                prior_sp = NULL,
                priors = NULL,
                constant_prior = TRUE,
                shift = 0.45, slope = 0.025,
                prior_sample_size = 5,
                min_prior_mean = 0.05,
                max_prior_mean = 0.95,
                summary_probs = c(0.025, 0.975),
                external_prevalence_data = NULL,
                prior_only = FALSE,
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
                                       thresholds = thresholds,
                                       prior_only = prior_only)
  if (is.null(priors)) {
    priors <- .get_prior_parameters(
      thresholds = thresholds,
      constant = constant,
      prior_p = prior_p,
      prior_se = prior_se,
      prior_sp = prior_sp,
      n_models_or_tests = n_models_or_tests,
      shift = shift,
      slope = slope,
      prior_sample_size = prior_sample_size,
      min_mean_se = min_prior_mean,
      max_mean_se = max_prior_mean,
      min_mean_sp = min_prior_mean,
      max_mean_sp = max_prior_mean
    )
  }
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

  if (is.null(external_prevalence_data)) {

    d_ext <- N_ext <- rep(0, n_thresholds)

  } else {

    if (any(external_prevalence_data < 0 | external_prevalence_data %% 1 != 0)) {
      stop("`external_prevalence_data` must be a vector with two non-negative integers")
    }
    if (external_prevalence_data[1] > external_prevalence_data[2]) {
      stop("`external_prevalence_data[1]` (cases) must be no greater than `external_prevalence_data[2]` (sample size)")
    }
    if (external_prevalence_data[2] < 1) {
      stop("`external_prevalence_data[2]` (sample size) must be at least 1.")
    }

    d_ext <- rep(external_prevalence_data[1], n_thresholds)
    N_ext <- rep(external_prevalence_data[2], n_thresholds)

  }
  other_models_indices <- lapply(
    1:n_models_or_tests,
    function(i)  (1:n_models_or_tests)[-i]
  )
  fit <- .dca_stan_list(
    n_thr = n_thresholds,
    n_models_or_tests = n_models_or_tests,
    N = rep(ifelse(prior_only, 0, N), n_thresholds),
    d = rep(ifelse(prior_only, 0, d), n_thresholds),
    tp = tp,
    tn = tn,
    thresholds = thresholds,
    N_ext = N_ext,
    d_ext = d_ext,
    prior_p1 = priors[['p1']],
    prior_p2 = priors[['p2']],
    prior_Se1 = priors[['Se1']],
    prior_Se2 = priors[['Se2']],
    prior_Sp1 = priors[['Sp1']],
    prior_Sp2 = priors[['Sp2']],
    other_models_indices = other_models_indices,
    refresh = refresh,
    ...
  )

  dca_summary <- .extract_dca_summary(fit = fit,
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
    output_data[['draws']] <- .extract_dca_draws(fit = fit,
                                                 model_or_test_names = model_or_test_names)
  }

  .output <- structure(output_data, class = "BayesDCAList")
  return(.output)
}

#' @title Get Threshold Performance Data for DCA list
#'
#' @param outcomes Integer vector (0 or 1) with binary outcomes.
#' @param predictions Numeric vector with predicted probabilities.
#' @param prior_only If TRUE, returns all zeros to made prior DCA.
#' @importFrom magrittr %>%
#' @keywords internal
.get_thr_data_list <- function(.data,
                               thresholds = seq(0.01, 0.5, 0.01),
                               prior_only = FALSE) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }
  outcomes <- .data[['outcomes']]
  model_or_test_names <- colnames(.data)[-1]
  N <- ifelse(prior_only, 0, length(outcomes))
  d <- ifelse(prior_only, 0, sum(outcomes))

  thr_data <- purrr::map(model_or_test_names, ~ {
    .predictions <- .data[[.x]]
    .thr_data <- tibble::tibble(
      model_or_test = .x,
      N = N, d = d, thresholds = thresholds
    ) %>%
      dplyr:::mutate(
        thr_perf = purrr::map(thresholds, function(.thr) {

          # # TODO: does the commented code below make sense?
          # if (.thr > 0.0) {
          #   # for binary tests, this gives the same tp/tn for all thresholds
          #   tp <- sum(.predictions[outcomes == 1] >= .thr)
          #   tn <- sum(.predictions[outcomes == 0] < .thr)
          # } else {
          #   tp <- sum(outcomes == 1)
          #   tn <- 0
          # }

          # for binary tests, this gives the same tp/tn for all thresholds
          # except thr = 0 in which case it gives `d`
          if (isTRUE(prior_only)) {
            tp <- 0
            tn <- 0
          } else {
            tp <- sum(.predictions[outcomes == 1] >= .thr)
            tn <- sum(.predictions[outcomes == 0] < .thr)
          }

          return(list(tp = tp, tn = tn))
        })
      ) %>%
      tidyr::unnest_wider(col = thr_perf) %>%
      dplyr::select(model_or_test, N, d, tp, tn, thresholds)
  }) %>%
    dplyr::bind_rows()

  return(thr_data)
}

#' @title Get summary from BayesDCA fit
#'
#' @param fit A stanfit object.
#' @param model_or_test_names Vector of names of models or binary tests under assessment.
#' @param summary_probs Numeric vector giving probabilities for credible interval.
#' @param thresholds Vector of thresholds for DCA.
#' @importFrom magrittr %>%
#' @keywords internal
.extract_dca_summary <- function(fit,
                                 model_or_test_names,
                                 summary_probs,
                                 thresholds) {

  # overall summary
  fit_summary <- rstan::summary(
    fit, probs = summary_probs
  )$summary %>%
    tibble::as_tibble(rownames = "par_name") %>%
    dplyr::select(-se_mean) %>%
    dplyr::rename(estimate := mean) %>%
    dplyr::mutate(
      threshold_ix = stringr::str_extract(par_name, "\\[\\d+") %>%
        stringr::str_remove(string = ., pattern = "\\[") %>%
        as.integer(),
      model_or_test_ix = stringr::str_extract(par_name, ",\\d+\\]") %>%
        stringr::str_remove_all(string = ., pattern = "\\]|,") %>%
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
      specificity = specificity
    ), class = "BayesDCASummary"
  )
  return(.summary)
}

#' @title Get posterior draws from DCA stanfit
#'
#' @param fit A stanfit object.
#' @param model_or_test_names Vector of names of models or binary tests under assessment.
#' @keywords internal
.extract_dca_draws <- function(fit,
                               model_or_test_names) {

  stan_draws <- rstan::extract(fit)

  .draws <- list(
    p = stan_draws$p,
    treat_all = stan_draws$treat_all,
    net_benefit = list(),
    delta_default = list(),
    prob_best = list(),
    Se = list(),
    Sp = list()
  )

  for (i in seq_along(model_or_test_names)) {
    .name <- model_or_test_names[i]
    .draws[['net_benefit']][[.name]] <- stan_draws$net_benefit[,,i]
    .draws[['delta_default']][[.name]] <- stan_draws$delta[,,i]
    .draws[['prob_best']][[.name]] <- stan_draws$prob_better_than_soc[,,i]
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

#' @title Get Threshold Performance Data
#'
#' @param outcomes Integer vector (0 or 1) with binary outcomes.
#' @param predictions Numeric vector with predicted probabilities.
#' @importFrom magrittr %>%
#' @keywords internal
get_thr_data <- function(outcomes,
                         predictions,
                         thresholds = seq(0.01, 0.5, 0.01)) {

  thr_data <- tibble::tibble(
    N = length(outcomes),
    d = sum(outcomes),
    thresholds = thresholds
  ) %>%
    dplyr:::mutate(
      thr_perf = purrr::map(thresholds, function(.thr) {
        tp <- sum(predictions[outcomes == 1] >= .thr)
        tn <- sum(predictions[outcomes == 0] < .thr)
        return(list(tp = tp, tn = tn))
      })
    ) %>%
    tidyr::unnest_wider(col = thr_perf) %>%
    dplyr::select(N, d, tp, tn, thresholds)
  return(thr_data)
}


