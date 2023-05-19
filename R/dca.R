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
#' @param summary_probs Probabilities used to compute credible
#' intervals (defaults to a 95% Cr.I.).
#' @param external_prevalence_data Vector with two positive
#' integers giving number of diseased and
#' non-diseased individuals, respectively, from external data
#' (e.g., if analyzing nested case-control data,
#' this is the number of cases and non-cases in the source
#' population).
#' @param refresh Control verbosity of `rstan::sampling`
#' (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter,
#' chains).
#' @return An object of class `BayesDCA`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot(fit)
dca <- function(.data,
                thresholds = seq(0, 0.5, length = 51),
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
                n_draws = 4000,
                ...) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }
  # avoid thresholds in {0, 1}
  thresholds <- thresholds %>%
    pmin(0.99) %>%
    pmax(1e-9) %>%
    unique()
  strategies <- colnames(.data)[-1]
  N <- nrow(.data)
  d <- sum(.data[["outcomes"]])
  n_thresholds <- length(thresholds)
  n_strategies <- length(strategies)
  threshold_data <- .get_thr_data_list(
    .data = .data,
    thresholds = thresholds,
    prior_only = prior_only
  )
  if (is.null(priors)) {
    priors <- .get_prior_parameters(
      thresholds = thresholds,
      constant = constant_prior,
      prior_p = prior_p,
      prior_se = prior_se,
      prior_sp = prior_sp,
      n_strategies = n_strategies,
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
    dplyr::select(thresholds, decision_strategy, tp) %>% # nolint
    tidyr::pivot_wider(
      names_from = decision_strategy,
      values_from = tp
    ) %>%
    dplyr::select(-thresholds) %>%
    as.matrix()

  tn <- threshold_data %>%
    dplyr::select(thresholds, decision_strategy, tn) %>% # nolint
    tidyr::pivot_wider(
      names_from = decision_strategy,
      values_from = tn
    ) %>%
    dplyr::select(-thresholds) %>%
    as.matrix()

  if (is.null(external_prevalence_data)) {
    d_ext <- N_ext <- 0
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

    d_ext <- external_prevalence_data[1]
    N_ext <- external_prevalence_data[2]
  }

  fit <- .dca_binary(
    n_thr = n_thresholds,
    strategies = strategies,
    N = ifelse(prior_only, 0, N),
    d = ifelse(prior_only, 0, d),
    tp = tp,
    tn = tn,
    thresholds = thresholds,
    N_ext = N_ext,
    d_ext = d_ext,
    prior_p1 = priors[["p1"]],
    prior_p2 = priors[["p2"]],
    prior_Se1 = priors[["Se1"]],
    prior_Se2 = priors[["Se2"]],
    prior_Sp1 = priors[["Sp1"]],
    prior_Sp2 = priors[["Sp2"]],
    n_draws = n_draws
  )

  dca_summary <- .extract_dca_summary(
    fit = fit,
    summary_probs = summary_probs,
    thresholds = thresholds,
    strategies = strategies
  )

  output_data <- list(
    fit = fit,
    summary = dca_summary,
    thresholds = thresholds,
    .data = .data,
    threshold_data = threshold_data,
    priors = priors,
    strategies = strategies
  )

  .output <- structure(output_data, class = "BayesDCA")
  return(.output)
}

#' Fit Bayesian Decision Curve Analysis using
#' Stan for list of models or binary tests
#'
#' @param n_thr Number of thresholds (int.).
#' @param n_strategies Number of models or binary tests (int.).
#' @param N Sample size (vector of integers of length `n_thr`).
#' @param d Diseased: number of diseased persons or
#' events (vector of integers of length `n_thr`).
#' @param tp True Positives: number of diseased persons correctly
#' identified as such by the diagnostic test of prediction
#' model (matrix of integers of size `n_thr` by `n_strategies`).
#' @param tn True Negatives: number of diseased persons correctly
#' identified as such by the diagnostic test of prediction
#' model (matrix of integers of size `n_thr` by
#' `n_strategies`).
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0.01, 0.5, 0.01)`).
#' @param N_ext,d_ext External sample size and number of
#' diseased individuals (or cases), respectively, used to
#' adjust prevalence.
#' @param prior_p,prior_se,prior_sp Prior parameters for
#' prevalence, sensitivity, and specificity (numeric matrices
#' of size `n_thr` by `n_strategies`).
#' @param refresh Control verbosity of
#' [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to
#' [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class
#' [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#' @keywords internal
.dca_binary <- function(n_thr,
                        strategies,
                        N,
                        d,
                        tp,
                        tn,
                        thresholds,
                        prior_p1, prior_p2,
                        prior_Se1, prior_Se2,
                        prior_Sp1, prior_Sp2,
                        other_models_indices,
                        N_ext = 0,
                        d_ext = 0,
                        n_draws = 4000) {
  thresholds <- pmin(thresholds, 0.999) # odds(1) = Inf
  n_strategies <- length(strategies)
  # get posterior distributions

  ## prevalence
  ### compute posterior parameters
  if (N_ext > 0) {
    stopifnot("You provided N_ext but d_ext is missing" = d_ext > 0)
    post_shape1_p <- d_ext + prior_p1
    post_shape2_p <- N_ext - d_ext + prior_p2
  } else {
    post_shape1_p <- d + prior_p1
    post_shape2_p <- N - d + prior_p2
  }
  ### sample from posterior
  p <- rbeta(n = n_draws, shape1 = post_shape1_p, shape2 = post_shape2_p)

  ## sensitivity, specificity, and all the rest
  ### initialize parameter variables
  post_shape1_Se <- post_shape2_Se <-
    post_shape1_Sp <- post_shape2_Sp <- matrix(
      nrow = n_thr, ncol = n_strategies,
      dimnames = list(
        NULL,
        strategies
      )
    )
  ### initialize distribution variables
  se <- sp <- net_benefit <- strategies %>%
    setNames(nm = .) %>%
    lapply(function(...) matrix(nrow = n_draws, ncol = n_thr))
  treat_all <- matrix(nrow = n_draws, ncol = n_thr)

  for (i in seq_len(n_thr)) {
    w_t <- thresholds[i] / (1 - thresholds[i])
    treat_all[, i] <- 1 * p - (1 - p) * (1 - 0) * w_t
    for (j in seq_len(n_strategies)) {
      # compute posterior parameters
      .m <- strategies[j]
      post_shape1_Se[i, j] <- tp[i, j] + prior_Se1[i, j]
      post_shape2_Se[i, j] <- d - tp[i, j] + prior_Se2[i, j]
      post_shape1_Sp[i, j] <- tn[i, j] + prior_Sp1[i, j]
      post_shape2_Sp[i, j] <- N - d - tn[i, j] + prior_Sp2[i, j]
      # sample from posteriors
      se[[.m]][, i] <- rbeta(
        n = n_draws,
        shape1 = post_shape1_Se[i, j],
        shape2 = post_shape2_Se[i, j]
      )
      sp[[.m]][, i] <- rbeta(
        n = n_draws,
        shape1 = post_shape1_Sp[i, j],
        shape2 = post_shape2_Sp[i, j]
      )

      # net benefit and other estimands
      net_benefit[[.m]][, i] <- se[[.m]][, i] * p - w_t * (1 - sp[[.m]][, i]) * (1 - p)
    }
  }

  # iterate again to look at "best competitor strategy"
  ## initialize comparison variables
  delta_useful <- delta_best <- strategies %>%
    setNames(nm = .) %>%
    lapply(function(...) matrix(nrow = n_draws, ncol = n_thr))
  p_useful <- p_best <- matrix(
    nrow = n_thr, ncol = n_strategies,
    dimnames = list(
      NULL,
      strategies
    )
  )
  for (i in seq_len(n_thr)) {
    nb_best_default_strategy <- pmax(0, treat_all[, i])
    for (j in seq_len(n_strategies)) {
      .m <- strategies[j]
      if (n_strategies == 1) {
        nb_best_competitor <- nb_best_default_strategy
      } else {
        nb_best_competitor <- matrixStats::rowMaxs(
          cbind(
            nb_best_default_strategy,
            sapply(net_benefit[-j], function(z) z[, i])
          )
        )
      }
      stopifnot(length(nb_best_competitor) == n_draws)
      delta_useful[[.m]][, i] <- net_benefit[[.m]][, i] - nb_best_default_strategy
      delta_best[[.m]][, i] <- net_benefit[[.m]][, i] - nb_best_competitor
      p_useful[i, j] <- mean(delta_useful[[.m]][, i] > 0)
      p_best[i, j] <- mean(delta_best[[.m]][, i] > 0)
    }
  }

  fit <- list(
    parameters = list(
      prevalence = list(shape1 = post_shape1_p, shape2 = post_shape2_p),
      sensitivity = list(shape1 = post_shape1_Se, shape2 = post_shape2_Se),
      specificity = list(shape1 = post_shape1_Se, shape2 = post_shape2_Se)
    ),
    distributions = list(
      prevalence = p,
      sensitivity = se,
      specificity = sp,
      net_benefit = net_benefit,
      treat_all = treat_all
    ),
    comparisons = list(
      delta_best = delta_best, p_best = p_best,
      delta_useful = delta_useful, p_useful = p_useful
    )
  )
  return(fit)
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
  outcomes <- .data[["outcomes"]]
  strategies <- colnames(.data)[-1]
  N <- ifelse(prior_only, 0, length(outcomes)) # nolint
  d <- ifelse(prior_only, 0, sum(outcomes)) # nolint

  thr_data <- purrr::map(
    strategies, ~ {
      .predictions <- .data[[.x]]
      .thr_data <- tibble::tibble(
        decision_strategy = .x,
        N = N, d = d, thresholds = thresholds
      ) %>%
        dplyr:::mutate(
          thr_perf = purrr::map(
            thresholds, function(.thr) {
              if (isTRUE(prior_only)) {
                tp <- 0
                tn <- 0
              } else {
                tp <- sum(.predictions[outcomes == 1] >= .thr)
                tn <- sum(.predictions[outcomes == 0] < .thr)
              }
              return(list(tp = tp, tn = tn))
            }
          )
        ) %>%
        tidyr::unnest_wider(col = thr_perf) %>%
        dplyr::select(decision_strategy, N, d, tp, tn, thresholds)
    }
  ) %>%
    dplyr::bind_rows()

  return(thr_data)
}


#' @title Get summary from BayesDCA fit
#'
#' @param fit A stanfit object.
#' @param strategies Vector of names of models or binary tests under assessment.
#' @param summary_probs Numeric vector giving probabilities for credible interval.
#' @param thresholds Vector of thresholds for DCA.
#' @importFrom magrittr %>%
#' @keywords internal
.extract_dca_summary <- function(fit, # nolint
                                 strategies,
                                 summary_probs,
                                 thresholds) {
  summ_labels <- paste0(
    round(summary_probs * 100, 1), "%"
  )
  .get_summ <- function(.x, .probs, .par_name, .summ_labels, .thrs, .decision_strategy_name) {
    if (is.vector(.x)) {
      .x <- as.matrix(.x, ncol = 1)
    }
    .qs <- matrix(
      matrixStats::colQuantiles(.x, prob = .probs),
      ncol = 2
    )
    tibble::tibble(
      par_name = .par_name,
      threshold = .thrs,
      decision_strategy_name = .decision_strategy_name,
      estimate = colMeans(.x),
      !!.summ_labels[1] := .qs[, 1],
      !!.summ_labels[2] := .qs[, 2]
    )
  }
  # prevalence
  prevalence_summary <- .get_summ(
    .x = fit$distributions$p,
    .probs = summary_probs,
    .par_name = "p",
    .summ_labels = summ_labels,
    .thrs = NA_real_,
    .decision_strategy_name = NA_character_
  )

  # treat all
  treat_all_summary <- .get_summ(
    .x = fit$distributions$treat_all,
    .probs = summary_probs,
    .par_name = "treat_all",
    .summ_labels = summ_labels,
    .thrs = thresholds,
    .decision_strategy_name = NA_character_
  )

  # se, sp, nb, deltas
  summary_by_model <- purrr::map_dfr(
    seq_along(strategies),
    function(j) {
      se <- .get_summ(
        .x = fit$distributions$sensitivity[[j]],
        .probs = summary_probs,
        .par_name = "Se",
        .summ_labels = summ_labels,
        .thrs = thresholds,
        .decision_strategy_name = strategies[j]
      )
      sp <- .get_summ(
        .x = fit$distributions$specificity[[j]],
        .probs = summary_probs,
        .par_name = "Sp",
        .summ_labels = summ_labels,
        .thrs = thresholds,
        .decision_strategy_name = strategies[j]
      )
      nb <- .get_summ(
        .x = fit$distributions$net_benefit[[j]],
        .probs = summary_probs,
        .par_name = "net_benefit",
        .summ_labels = summ_labels,
        .thrs = thresholds,
        .decision_strategy_name = strategies[j]
      )
      d_best <- .get_summ(
        .x = fit$comparisons$delta_best[[j]],
        .probs = summary_probs,
        .par_name = "delta_best",
        .summ_labels = summ_labels,
        .thrs = thresholds,
        .decision_strategy_name = strategies[j]
      )
      d_useful <- .get_summ(
        .x = fit$comparisons$delta_best[[j]],
        .probs = summary_probs,
        .par_name = "delta_useful",
        .summ_labels = summ_labels,
        .thrs = thresholds,
        .decision_strategy_name = strategies[j]
      )
      x <- dplyr::bind_rows(
        nb, d_best, d_useful, se, sp
      )
    }
  ) %>%
    dplyr::arrange(par_name, threshold, decision_strategy_name)

  # sensitivity
  sensitivity <- summary_by_model %>%
    dplyr::filter(
      stringr::str_detect(par_name, "Se") # nolint
    )

  # specificity
  specificity <- summary_by_model %>%
    dplyr::filter(
      stringr::str_detect(par_name, "Sp") # nolint
    )

  # net benefit
  net_benefit <- summary_by_model %>%
    dplyr::filter(stringr::str_detect(par_name, "net_benefit")) # nolint

  # deltas
  delta_best <- summary_by_model %>%
    dplyr::filter(stringr::str_detect(par_name, "delta_best")) # nolint
  delta_useful <- summary_by_model %>%
    dplyr::filter(stringr::str_detect(par_name, "delta_useful")) # nolint


  .summary <- structure(
    list(
      net_benefit = net_benefit,
      treat_all = treat_all_summary,
      prevalence = prevalence_summary,
      sensitivity = sensitivity,
      specificity = specificity,
      delta_best = delta_best,
      delta_useful = delta_useful
    ),
    class = "BayesDCASummary"
  )
  return(.summary)
}


#' @title Print BayesDCA
#'
#' @param obj BayesDCA object
#' @export
print.BayesDCA <- function(obj, ...) {
  .data <- paste0(
    "N=", unique(obj$threshold_data$N), "\tD=", unique(obj$threshold_data$d)
  )
  cat(
    paste0(
      c(
        "BayesDCA\n",
        paste0("Number of thresholds: ", length(obj$thresholds)),
        "\nRaw data:", .data,
        "Models or tests: ",
        paste0(obj$strategies, collapse = ", ")
      ),
      collapse = "\n"
    )
  )
}

#' @title Get delta plot data for binary bayesDCA
#'
#' @param obj BayesDCA object
#' @param strategies Character vector with subset of decision strategies.
#' @param type One of "best", "useful", or "pairwise".
#' @importFrom magrittr %>%
#' @keywords internal
#' @return A ggplot object.
get_superiority_prob_plot_data_binary <- function(obj, # nolint: object_length_linter.
                                                  min_diff = 0,
                                                  strategies = NULL,
                                                  type = c("best", "useful", "pairwise"),
                                                  labels = NULL) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two strategies to plot pairwise comparison" = length(strategies) == 2 # nolint
    )
  }

  strategies <- validate_strategies(
    obj = obj, strategies = strategies
  )

  if (type == "pairwise") {
    nb1 <- obj$fit$distributions$net_benefit[[strategies[1]]]
    nb2 <- obj$fit$distributions$net_benefit[[strategies[2]]]
    .prob <- colMeans(nb1 - nb2 > min_diff)
    df <- tibble::tibble(
      prob = .prob,
      threshold = obj$thresholds
    )
  } else {
    df <- lapply(
      seq_along(strategies),
      function(j) {
        .m <- strategies[j]
        if (type == "best") {
          .prob <- obj$fit$comparisons$p_best[, j]
        } else {
          .prob <- obj$fit$comparisons$p_useful[, j]
        }
        tibble::tibble(
          prob = .prob,
          threshold = obj$thresholds,
          decision_strategy = .m,
          label = labels[.m]
        )
      }
    ) %>%
      dplyr::bind_rows()
  }
  return(df)
}

#' @title Get delta plot data for binary bayesDCA
#'
#' @param obj BayesDCA object
#' @param strategies Character vector with subset of decision strategies.
#' @param type One of "best", "useful", or "pairwise".
#' @importFrom magrittr %>%
#' @keywords internal
#' @return A ggplot object.
get_delta_plot_data_binary <- function(obj,
                                       strategies = NULL,
                                       type = c("best", "useful", "pairwise"),
                                       labels = NULL) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two strategies to plot pairwise comparison" = length(strategies) == 2 # nolint
    )
  }

  strategies <- validate_strategies(
    obj = obj, strategies = strategies
  )

  if (type == "pairwise") {
    nb1 <- obj$fit$distributions$net_benefit[[strategies[1]]]
    nb2 <- obj$fit$distributions$net_benefit[[strategies[2]]]
    .delta <- nb1 - nb2
    qs <- matrixStats::colQuantiles(.delta, probs = c(.025, .975))
    df <- tibble::tibble(
      estimate = colMeans(.delta),
      `2.5%` = qs[, "2.5%"],
      `97.5%` = qs[, "97.5%"],
      threshold = obj$thresholds,
    )
  } else {
    df <- lapply(
      seq_along(strategies),
      function(i) {
        .m <- strategies[i]
        if (type == "best") {
          .delta <- obj$fit$comparisons$delta_best[[.m]]
        } else {
          .delta <- obj$fit$comparisons$delta_useful[[.m]]
        }
        qs <- matrixStats::colQuantiles(.delta, probs = c(.025, .975)) # nolint
        tibble::tibble(
          estimate = colMeans(.delta),
          `2.5%` = qs[, "2.5%"],
          `97.5%` = qs[, "97.5%"],
          threshold = obj$thresholds,
          decision_strategy = .m,
          label = labels[.m]
        )
      }
    ) %>%
      dplyr::bind_rows()
  }
  return(df)
}
