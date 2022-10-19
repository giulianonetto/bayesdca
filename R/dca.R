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
                summary_probs = c(0.025, 0.975),
                external_prevalence_data = NULL,
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

  fit <- .dca_stan_list(
    n_thr = n_thresholds,
    n_models_or_tests = n_models_or_tests,
    N = rep(N, n_thresholds),
    d = rep(d, n_thresholds),
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
#' @importFrom magrittr %>%
.get_thr_data_list <- function(.data,
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
.get_prior_parameters <- function(n_thresholds,
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
.extract_dca_summary <- function(fit,
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
      specificity = specificity
    ), class = "BayesDCASummary"
  )
  return(.summary)
}

#' @title Get posterior draws from DCA stanfit
#'
#' @param fit A stanfit object.
#' @param model_or_test_names Vector of names of models or binary tests under assessment.
.extract_dca_draws <- function(fit,
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
    .draws[['delta_default']][[.name]] <- stan_draws$delta[,,i]
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

#' @title Plot BayesDCAList
#'
#' @param obj BayesDCAList object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @param models_or_tests Character vector with models or tests to compare. If null, compares either first two in `obj$model_or_test_names` or the first one against Treat all/none (if only one available).
#' @importFrom magrittr %>%
#' @export
plot.BayesDCAList <- function(obj,
                              models_or_tests = NULL,
                              colors = NULL,
                              labels = NULL,
                              raw_values = NULL,
                              raw_values_label = "Biomarker threshold", ...) {


  if (!is.null(models_or_tests)) {

    stopifnot(
      "Provided `models_or_tests` are not available" = all(
        models_or_tests %in% obj$model_or_test_names
      )
    )

    net_benefit_data <- obj$summary$net_benefit %>%
      dplyr::filter(model_or_test_name %in% models_or_tests)

  } else {

    net_benefit_data <- obj$summary$net_benefit

  }

  colors_and_labels <- get_colors_and_labels(obj = obj,
                                             colors = colors,
                                             labels = labels,
                                             models_or_tests = models_or_tests)

  .ymin <- ifelse(
    max(obj$summary$treat_all$estimate) > 0.02,
    -0.02,
    -max(obj$summary$treat_all$estimate)
  )

  .p <- ggplot2::ggplot() +
    # set x axis
    ggplot2::aes(x = threshold) +
    # add color/fill/label scheme
    colors_and_labels +
    # add treat all curve
    ggplot2::geom_ribbon(
      data = obj$summary$treat_all,
      ggplot2::aes(ymax = `97.5%`, ymin = `2.5%`,
                   fill = "Treat all"),
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = obj$summary$treat_all,
      ggplot2::aes(y = estimate, color = "Treat all")
    ) +
    # add net benefit curves
    ggplot2::geom_ribbon(
      data = net_benefit_data,
      ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`,
                   fill = model_or_test_name),
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = net_benefit_data,
      ggplot2::aes(y = estimate,
                   color = model_or_test_name)
    ) +
    # add treat none curve
    ggplot2::geom_hline(
      ggplot2::aes(color = "Treat none", yintercept = 0),
      linetype = 'longdash', lwd = 0.8
    ) +
    # make it pretty
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::coord_cartesian(ylim = c(.ymin, NA)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(x = "Decision threshold", y = "Net Benefit",
                  color = NULL) +
    ggplot2::guides(fill = 'none')

  if (!is.null(raw_values)) {

    stopifnot(is.data.frame(raw_values))

    raw_values_p <- ggplot2::ggplot(
      raw_values,
      ggplot2::aes(x = thresholds*100)
    ) +
      ggplot2::scale_x_continuous(
        labels = raw_values$values,
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        x = raw_values_label
      )
    .p <- (.p / raw_values_p) +
      patchwork::plot_layout(
        heights = c(0.9, 0.01)
      )
  }

  return(.p)
}

#' @title Plot BayesDCAList comparison
#'
#' @param obj BayesDCAList object
#' @param models_or_tests Character vector with models or tests to compare. If null, compares either first two in `obj$model_or_test_names` or the first one against Treat all/none (if only one available).
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @param plot_list If TRUE, returns a list of separate ggplot objects.
#' @importFrom magrittr %>%
#' @import patchwork
#' @export
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' compare_dca(fit)
#' @return A patchwork/ggplot object or a list of ggplot objects.
compare_dca <- function(obj, models_or_tests = NULL, colors = NULL, labels = NULL,
                        plot_list = FALSE, .evpi = FALSE, ...) {

  if (is.null(models_or_tests)) {
    models_or_tests <- as.vector(na.omit(obj$model_or_test_names[1:2]))
  } else {
    stopifnot(
      "Provided `models_or_tests` are not available" = all(
        models_or_tests %in% obj$model_or_test_names
      )
    )
  }

  stopifnot(length(models_or_tests) > 0 & length(models_or_tests) < 3)

  p1 <- plot(obj, colors = colors, labels = labels,
             models_or_tests = models_or_tests)
  p2 <- plot_delta_nb(obj = obj,
                      models_or_tests = models_or_tests,
                      labels = labels)
  p3 <- plot_prob_better(obj = obj,
                         models_or_tests = models_or_tests,
                         labels = labels)

  if (isTRUE(plot_list)) {
    .plot_list <- list(dca = p1,
                       delta = p2,
                       prob_better = p3)
    if (isTRUE(.evpi)) {
      .plot_list[['evpi']] <- plot_evpi(obj = obj,
                                        models_or_tests = models_or_tests,
                                        labels = labels)

    }

    return(.plot_list)

  } else {
    if (isTRUE(.evpi)) {
      p4 <- plot_evpi(obj = obj,
                      models_or_tests = models_or_tests,
                      labels = labels)
      .p <- (p1|p3)/(p2|p4) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = 'bottom')
    } else {
      .p <- p1/(p2 | p3)
    }
    return(.p)
  }
}

#' @title Plot BayesDCAList delta
#'
#' @param obj BayesDCAList object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot_delta_nb(fit)
#' @return A ggplot object.
plot_delta_nb <- function(obj, models_or_tests, labels = NULL) {

  if (is.null(obj$draws)) {
    stop("Missing draws in BayesDCAList object. Run dca() with `keep_draws = TRUE`")
  }

  # build labels for plot subtitle
  plot_labels <- vector("character", length = 2L)

  if (models_or_tests[1] %in% names(labels)) {
    plot_labels[1] <- labels[models_or_tests[1]]
  } else {
    plot_labels[1] <- models_or_tests[1]
  }

  if (length(models_or_tests) > 1) {
    if (models_or_tests[2] %in% names(labels)) {
      plot_labels[2] <- labels[models_or_tests[2]]
    } else {
      plot_labels[2] <- models_or_tests[2]
    }
  }


  # retrieve posterior draws
  if (length(models_or_tests) == 1) {
    # get delta posterior draws
    .delta <- obj$draws$delta_default[[models_or_tests]]
    # get subtitles
    .subtitle <- paste0(plot_labels, ' \u2212 Treat all or none')
  } else {
    # get delta posterior draws
    nb1 <- obj$draws$net_benefit[[models_or_tests[1]]]
    nb2 <- obj$draws$net_benefit[[models_or_tests[2]]]
    .delta <- nb1 - nb2
    # get subtitles
    .subtitle <- paste0(plot_labels[1],
                        ' \u2212 ',
                        plot_labels[2])
  }

  q <- matrixStats::colQuantiles(.delta, probs = c(.025, .5, .975))

  .plot <- tibble::tibble(
    threshold = obj$thresholds,
    estimate = q[,'50%'],
    `2.5%` = q[,'2.5%'],
    `97.5%` = q[,'97.5%']
  ) %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = threshold, y = estimate,
                 ymin = `2.5%`, ymax = `97.5%`) +
    ggplot2::geom_ribbon(alpha = 0.4) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 'longdash',
      color = 'gray30', lwd = 0.8
    ) +
    ggplot2::labs(
      x = "Decision threshold",
      y = '\u0394 NB',
      subtitle = .subtitle
    )
  return(.plot)
}

#' @title Plot probability of better net benefit between two models or tests
#' @param obj BayesDCAList object
#' @param models_or_tests Character vector with models or tests to compare. If null, compares either first two in `obj$model_or_test_names` or the first one against Treat all/none (if only one available).
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot_prob_better(fit)

#' @return A ggplot object.
plot_prob_better <- function(obj, models_or_tests = NULL, labels = NULL) {
  if (is.null(obj$draws)) {
    stop("Missing draws in BayesDCAList object. Run dca() with `keep_draws = TRUE`")
  }

  if (is.null(models_or_tests)) {
    models_or_tests <- as.vector(na.omit(obj$model_or_test_names[1:2]))
  } else {
    stopifnot(
      "Provided `models_or_tests` are not available" = all(
        models_or_tests %in% obj$model_or_test_names
      )
    )
  }

  # build labels for plot subtitle
  plot_labels <- vector("character", length = 2L)

  if (models_or_tests[1] %in% names(labels)) {
    plot_labels[1] <- labels[models_or_tests[1]]
  } else {
    plot_labels[1] <- models_or_tests[1]
  }

  if (length(models_or_tests) > 1) {
    if (models_or_tests[2] %in% names(labels)) {
      plot_labels[2] <- labels[models_or_tests[2]]
    } else {
      plot_labels[2] <- models_or_tests[2]
    }
  }

  # retrieve posterior draws
  if (length(models_or_tests) == 1) {
    # get delta posterior draws
    .delta <- obj$draws$delta_default[[models_or_tests]]
    # get subtitles
    .subtitle <- paste0("P(", plot_labels, ' \u003E  Treat all or none)')
  } else {
    # get delta posterior draws
    nb1 <- obj$draws$net_benefit[[models_or_tests[1]]]
    nb2 <- obj$draws$net_benefit[[models_or_tests[2]]]
    .delta <- nb1 - nb2
    # get ubtitles
    .subtitle <- paste0("P(",
                        plot_labels[1],
                        ' \u003E ',
                        plot_labels[2],
                        ")")
  }

  .plot <- tibble::tibble(
    threshold = obj$thresholds,
    estimate = colMeans(.delta > 0)
  ) %>%
    ggplot2::ggplot(ggplot2::aes(x = threshold, y = estimate)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(1),
      breaks = scales::pretty_breaks(10),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Decision threshold", y = NULL,
      subtitle = .subtitle
    )

  return(.plot)
}

#' @title Get Threshold Performance Data
#'
#' @param outcomes Integer vector (0 or 1) with binary outcomes.
#' @param predictions Numeric vector with predicted probabilities.
#' @importFrom magrittr %>%
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


