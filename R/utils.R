#' @title Compute Expected Value of Perfect Information (EVPI)
#'
#' @param obj BayesDCAList or BayesDCASurv object
#' @param ... matrices of dimension n_draws * n_thresholds
#' containing posterior net benefit
#' @importFrom magrittr %>%
#' @keywords internal
evpi <- function(thresholds, ...) {
  .dots <- list(...)
  .evpi <- numeric(length = length(thresholds))
  for (i in seq_along(thresholds)) {
    .dots_i <- lapply(.dots, function(.d) .d[, i])
    mean_nbs <- unlist(lapply(.dots_i, function(nb_draws) mean(nb_draws)))
    max_nb_draws <- matrixStats::rowMaxs(
      cbind(0, do.call(cbind, .dots_i))
    )
    ENB_perfect <- mean(max_nb_draws) # nolint
    ENB_current <- max(0, mean_nbs) # nolint
    .evpi[i] <- ENB_perfect - ENB_current
  }
  return(.evpi)
}

#' @title Minimal events per interval
#' @keywords internal
min_events_per_interval <- function() {
  return(5)
}

#' @title Get cutpoints for survival estimation
#' @param .prediction_time time point at which event is predicted to happen
#' @param .event_times times of observed events (non-censored)
#' @param .base_cutpoints vector of cutpoints to start with
#' @keywords internal
get_cutpoints <- function(.prediction_time,
                          .event_times,
                          .base_cutpoints = c(0.25, 0.5, 0.75, 1)) {
  stopifnot("All event times must be positive." = all(.event_times > 0))
  min_events <- min_events_per_interval()
  .base_cutpoints <- .base_cutpoints[.base_cutpoints > 0] * .prediction_time
  events_above_cutpoint <- sapply(
    .base_cutpoints, function(cutpoint) sum(.event_times > cutpoint)
  )
  .base_cutpoints <- .base_cutpoints[events_above_cutpoint >= min_events]
  # only keep cutpoints that correspond to
  # intervals with at least `min_events` events
  .previous <- new_cutpoints <- 0
  for (i in seq_along(.base_cutpoints)) {
    .current <- .base_cutpoints[i]
    events <- sum(.event_times > .previous & .event_times <= .current)
    if (events >= min_events) {
      new_cutpoints <- c(new_cutpoints, .current)
      .previous <- .current
    }
  }

  return(new_cutpoints)
}

#' @title Get events per intervals defined by set of cutpoints
#' @param .cutpoints cutpoints defining intervals
#' @param .event_times times of observed events (non-censored)
#' @keywords internal
get_events_per_interval <- function(.cutpoints, .event_times) {
  table(cut(.event_times, c(.cutpoints, Inf)))
}


#' @title Get colors and labels for BayesDCA plots
#'
#' @param obj BayesDCAList object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @keywords internal
get_colors_and_labels <- function(obj,
                                  models_or_tests = NULL,
                                  colors = NULL, labels = NULL,
                                  all_or_none = TRUE) {
  # decide which models/tests to include
  if (is.null(models_or_tests)) {
    model_or_test_names <- obj$model_or_test_names
  } else {
    stopifnot(
      any(models_or_tests %in% obj$model_or_test_names)
    )
    model_or_test_names <- models_or_tests[
      models_or_tests %in% obj$model_or_test_names
    ]
  }
  # pick color palette for ggplot
  if (isTRUE(all_or_none)) {
    color_values <- c(
      "Treat all" = "black", "Treat none" = "gray40"
    )
  } else {
    color_values <- character()
  }

  n_colors <- length(model_or_test_names)
  if (n_colors < 9) {
    palette <- RColorBrewer:::brewer.pal(max(c(n_colors, 3)), "Dark2")
  } else {
    palette <- grDevices::colorRampPalette(
      RColorBrewer:::brewer.pal(n_colors, "Set2")
    )(n_colors)
  }
  # set actual color values to use in scale_color_manual
  for (i in seq_len(n_colors)) {
    model_or_test <- model_or_test_names[i]
    if (!is.null(colors) && model_or_test %in% names(colors)) {
      color_values[model_or_test] <- colors[[model_or_test]]
    } else {
      color_values[model_or_test] <- palette[i]
    }
  }
  # define color and label scales
  if (is.null(labels)) {
    colors_and_labels <- list(
      ggplot2::scale_color_manual(
        values = color_values
      ),
      ggplot2::scale_fill_manual(
        values = color_values
      )
    )
  } else {
    colors_and_labels <- list(
      ggplot2::scale_color_manual(
        labels = labels, values = color_values
      ),
      ggplot2::scale_fill_manual(
        labels = labels, values = color_values
      )
    )
  }

  return(colors_and_labels)
}

#' @title Get time exposed within each interval for a given prediction time
#'
#' @param .prediction_time Time for event prediction
#' @param .cutpoints Cutpoints for constant hazard interval
#' @keywords internal
get_survival_time_exposed <- function(.prediction_time, .cutpoints) {
  time_exposed <- numeric(length = length(.cutpoints))
  for (i in seq_len(length(.cutpoints))) {
    .lower <- .cutpoints[i]
    .upper <- c(.cutpoints, Inf)[i + 1]
    if (.prediction_time >= .upper) {
      # if prediction time > upper bound, use interval length
      time_exposed[i] <- .upper - .lower
    } else if (.prediction_time > .lower) {
      # if prediction time > lower bound,
      # use distance between pred time and lower bound
      time_exposed[i] <- .prediction_time - .lower
    } else {
      # otherwise, prediction time is before interval, exposure time is zero
      time_exposed[i] <- 0
    }
  }
  return(time_exposed)
}

#' @title Get death pseudo-counts
#'
#' @param .prediction_data Contains prognostic model predictions
#' @param .surv_data Contains observed time to event
#' (`.time` column) and event indicator (`.status` column)
#' @param .cutpoints Survival cutpoints.
#' @param .models_or_tests Names for models or tests.
#' @param .thresholds DCA thresholds.
#' @param .prior_scaling_factor Prior
#' alpha = (0.69/median_surv)*.prior_scaling_factor,
#' Prior beta = .prior_scaling_factor.
#' @importFrom magrittr %>%
#' @importFrom survival Surv
#' @keywords internal
get_death_pseudo_counts <- function(.prediction_data, # nolint
                                    .surv_data,
                                    .cutpoints,
                                    .models_or_tests,
                                    .thresholds) {
  n_cuts <- length(.cutpoints)
  n_thr <- length(.thresholds)
  death_pseudo_counts <- vector("list", length = length(.models_or_tests))
  total_exposure_times <- vector("list", length = length(.models_or_tests))

  for (i in seq_along(.models_or_tests)) {
    .model <- .models_or_tests[i]
    death_pseudo_counts[[i]] <- matrix(nrow = n_cuts, ncol = n_thr)
    total_exposure_times[[i]] <- matrix(nrow = n_cuts, ncol = n_thr)
    for (j in seq_along(.thresholds)) {
      .thr <- .thresholds[j]
      if (all(is.na(.prediction_data))) {
        .positive_prediction <- 1 >= .thr # for S0, all positives
      } else {
        .predictions <- .prediction_data[[.model]]
        .positive_prediction <- .predictions >= .thr
      }

      .d <- .surv_data[.positive_prediction, ]
      stopifnot("no positive prediction not implemented" = nrow(.d) > 0)
      .d$patient_id <- 1:nrow(.d) # nolint
      .d_split <- survival::survSplit(
        Surv(.time, .status) ~ 1, # nolint
        data = .d,
        cut = .cutpoints,
        subset = 1:nrow(.d), # nolint
        id = "patient_id",
        start = "tstart",
        end = "tstop"
      ) %>%
        dplyr::group_by(
          interval_id = paste0(
            "interval_", as.numeric(factor(tstart)) # nolint
          )
        ) %>%
        dplyr::summarise(
          total_events = sum(.status),
          total_exposure_time = sum(tstop - tstart) # nolint
        )
      death_pseudo_counts[[i]][, j] <- .d_split$total_events
      total_exposure_times[[i]][, j] <- .d_split$total_exposure_time
    }
  }
  .output <- list(
    death_pseudo_counts = death_pseudo_counts,
    total_exposure_times = total_exposure_times
  )
  return(.output)
}

#' @title Get posterior parameters for survival model
#'
#' @param .prediction_data Contains prognostic model predictions
#' @param .surv_data Contains observed time to event
#' (`.time` column) and event indicator (`.status` column)
#' @param .cutpoints Survival cutpoints.
#' @param .models_or_tests Names for models or tests.
#' @param .thresholds DCA thresholds.
#' @param .prior_scaling_factor Prior
#' @param .prior_only If TRUE, will sample from prior only.
#' alpha = (0.69/median_surv)*.prior_scaling_factor,
#' Prior beta = .prior_scaling_factor.
#' @importFrom magrittr %>%
#' @importFrom survival Surv
#' @keywords internal
get_survival_posterior_parameters <- function(.prediction_data, # nolint
                                              .surv_data,
                                              .cutpoints,
                                              .models_or_tests,
                                              .thresholds,
                                              .prior_scaling_factor,
                                              .prediction_time,
                                              .prior_anchor = c("median", "prediction_time"), # nolint
                                              .prior_only = FALSE,
                                              .prior_means = NULL) {
  .prior_anchor <- match.arg(.prior_anchor)
  n_cuts <- length(.cutpoints)
  n_thr <- length(.thresholds)
  initialize_pars <- function() {
    lapply(
      seq_along(.models_or_tests),
      function(...) {
        matrix(
          nrow = n_cuts,
          ncol = n_thr
        )
      }
    )
  }

  all_posterior_alphas <- initialize_pars()
  all_posterior_betas <- initialize_pars()

  for (i in seq_along(.models_or_tests)) {
    .model <- .models_or_tests[i]
    w_inv <- .prior_scaling_factor # might be updated for some thresholds
    empty_thresholds <- 0
    for (j in seq_along(.thresholds)) {
      .thr <- .thresholds[j]
      if (all(is.na(.prediction_data))) {
        .positive_prediction <- 1 >= .thr # for S0, all positives
      } else {
        .predictions <- .prediction_data[[.model]]
        .positive_prediction <- .predictions >= .thr
      }

      .d <- .surv_data[.positive_prediction, ]

      if (is.null(.prior_means)) { # use .prior_anchor for prior mean
        .n_events <- sum(.d$.status == 1L)
        if (.n_events > 0) {
          if (.prior_anchor == "median") {
            .median_surv <- survival:::median.Surv(
              Surv(.d$.time, .d$.status)
            )$quantile
            .prior_mean <- -log(0.5) / .median_surv
          } else { # prior anchor = "prediction_time"
            max_observed_time <- max(.d$.time)
            if (.prediction_time > max_observed_time) {
              msg <- paste0(
                "Prediction time (",
                .prediction_time,
                ") is greater than the largest observed time (",
                max_observed_time,
                ") for threshold ", .thr,
                ".\n Either use decision thresholds lower than ",
                .thr,
                " or set .prior_anchor = 'median'."
              )
              stop(msg)
            }
            obs_surv <- survival:::summary.survfit(
              survival:::survfit(Surv(.d$.time, .d$.status) ~ 1),
              time = .prediction_time,
              extend = TRUE
            )[["surv"]]
            .prior_mean <- -log(obs_surv) / .prediction_time
          }
        } else {
          msg <- paste0(
            "Zero events among positive patients for model '",
            .model, "' and threshold ", .thr,
            ". Using median survival as prior anchor."
          )
          message(cli::col_red(msg))
          empty_thresholds <- empty_thresholds + 1
          # positive patients with previous threshold
          thr_ix <- j - empty_thresholds
          .d_previous <- .surv_data[.predictions >= .thresholds[thr_ix], ]
          .median_surv <- survival:::median.Surv(
            Surv(.d_previous$.time, .d_previous$.status)
          )$quantile
          .prior_mean <- -log(0.5) / .median_surv
          w_inv <- w_inv / 2
        }

        if (is.na(.prior_mean)) {
          msg <- "Failed to compute prior mean. Setting prior mean to 1."
          message(cli::col_red(msg))
          .prior_mean <- 1
        }
      } else {
        .prior_mean <- .prior_means[j]
      }
      .prior_alpha <- .prior_mean * w_inv
      .prior_beta <- w_inv

      if (nrow(.d) > 0 && isFALSE(.prior_only)) {
        .d$patient_id <- 1:nrow(.d) # nolint
        .d_split <- survival::survSplit(
          Surv(.time, .status) ~ 1, # nolint
          data = .d,
          cut = .cutpoints,
          subset = 1:nrow(.d), # nolint
          id = "patient_id",
          start = "tstart",
          end = "tstop"
        ) %>%
          # TODO: might want to double check the calculation below
          dplyr::group_by(
            interval_id = paste0(
              "interval_", as.numeric(factor(tstart)) # nolint
            )
          ) %>%
          dplyr::summarise(
            total_events = sum(.status),
            total_exposure_time = sum(tstop - tstart) # nolint
          )

        n_empty_intervals <- n_cuts - nrow(.d_split)
        if (n_empty_intervals > 0) {
          msg <- paste0(
            "Got ", n_empty_intervals, " empty interval (s) for threshold ",
            .thresholds[j], " in model '", .model, "'"
          )
          message(cli::col_cyan(msg))
          for (k in 1:n_empty_intervals) {
            .d_split <- .d_split %>%
              dplyr::add_row(
                interval_id = paste0("interval_", nrow(.d_split) + k),
                total_events = 0,
                total_exposure_time = 0
              )
          }
        }

        all_posterior_alphas[[i]][, j] <- .d_split$total_events + .prior_alpha
        all_posterior_betas[[i]][, j] <- .d_split$total_exposure_time + .prior_beta # nolint
      } else {
        all_posterior_alphas[[i]][, j] <- 0 + .prior_alpha
        all_posterior_betas[[i]][, j] <- 0 + .prior_beta
      }
    }
    if (empty_thresholds > 0) {
      msg <- paste0(
        "Observed ", empty_thresholds,
        " thresholds with zero events among positive patients for model '",
        .model, "'"
      )
      message(cli::col_br_red(msg))
    }
  }

  if (length(.thresholds) == 1) {
    all_posterior_alphas <- unlist(all_posterior_alphas)
    all_posterior_betas <- unlist(all_posterior_betas)
  }

  .posterior_pars <- list(
    .alpha = all_posterior_alphas,
    .beta = all_posterior_betas
  )
  return(.posterior_pars)
}


#' @title Get posterior parameters for positivity probability
#'
#' @param .prediction_data Contains prognostic model predictions
#' @param .thresholds DCA thresholds.
#' @param .prior_shape1 Shape 1 for beta prior.
#' @param .prior_shape2 Shape 2 for beta prior.
#' @keywords internal
get_positivity_posterior_parameters <- function(.prediction_data, # nolint
                                                .thresholds,
                                                .prior_shape1 = 1,
                                                .prior_shape2 = 1,
                                                .prior_only = FALSE) {
  N <- nrow(.prediction_data) # nolint
  n_models <- ncol(.prediction_data)
  .models_or_tests <- colnames(.prediction_data)
  n_thresholds <- length(.thresholds)

  all_posterior_shape1 <- matrix(
    nrow = n_models,
    ncol = n_thresholds
  )
  all_posterior_shape2 <- matrix(
    nrow = n_models,
    ncol = n_thresholds
  )

  for (i in seq_along(.models_or_tests)) {
    .model <- .models_or_tests[i]
    for (j in seq_along(.thresholds)) {
      if (isFALSE(.prior_only)) {
        .thr <- .thresholds[j]
        .predictions <- .prediction_data[[.model]]
        .positive_prediction <- .predictions >= .thr
        total_positives <- sum(.positive_prediction)
        all_posterior_shape1[i, j] <- total_positives + .prior_shape1
        all_posterior_shape2[i, j] <- N - total_positives + .prior_shape2
      } else {
        all_posterior_shape1[i, j] <- .prior_shape1
        all_posterior_shape2[i, j] <- .prior_shape2
      }
    }
  }

  .posterior_pars <- list(
    .shape1 = all_posterior_shape1,
    .shape2 = all_posterior_shape2
  )
  return(.posterior_pars)
}

#' Get threshold-specific mean prior sensitivity
#' @param thresholds Vector of decision thresholds.
#' @param shift Scalar controlling height of prior Sensitivity curve
#' @param slope Scalar controlling shape of prior Sensitivity curve
#' @param .min,.max Minimum and maximum prior mean
#' @importFrom magrittr %>%
#' @keywords internal
get_prior_se_mu <- function(thresholds, shift = 0.45,
                            slope = 0.025, .min = 0.1, .max = 0.9) {
  x <- plogis(shift + slope * qlogis(1 - thresholds)^3)
  x %>%
    # prior mean cannot be exactly 0 or 1
    pmax(.min) %>%
    pmin(.max)
}

#' Get threshold-specific mean prior specificity
#' @param thresholds Vector of decision thresholds.
#' @param shift Scalar controlling height of prior Specificity curve
#' @param slope Scalar controlling shape of prior Specificity curve
#' @param .min,.max Minimum and maximum prior mean
#' @importFrom magrittr %>%
#' @keywords internal
get_prior_sp_mu <- function(thresholds, shift = 0.45,
                            slope = 0.025, .min = 0.1, .max = 0.9) {
  x <- plogis(shift + slope * qlogis(thresholds)^3)
  x %>%
    # prior mean cannot be exactly 0 or 1
    pmax(.min) %>%
    pmin(.max)
}

#' Get threshold-specific prior sample size or strength
#' @param min_prior_sample_size Minimum prior sample size or strength.
#' @param max_prior_sample_size Maximum prior sample size or strength.
#' @param slope_prior_sample_size Rate of change in
#' prior sample size or strength.
#' @param thresholds Decision thresholds
#' @keywords internal
get_prior_sample_size <- function(thresholds,
                                  min_prior_sample_size = 20,
                                  max_prior_sample_size = 100,
                                  slope_prior_sample_size = 300) {
  x <- (
    max_prior_sample_size
    - slope_prior_sample_size * thresholds
      + slope_prior_sample_size * thresholds^2
  )
  x %>%
    pmax(min_prior_sample_size) %>%
    pmin(max_prior_sample_size)
}

#' @title Get priors for Bayesian DCA
#'
#' @param thresholds Vector of decision thresholds.
#' @param shift Scalar controlling height of prior
#' Specificity curve. Only used if `constant=FALSE`.
#' @param slope Scalar controlling shape of prior
#' Specificity curve. Only used if `constant=FALSE`.
#' @param min_mean_se,min_mean_sp,max_mean_se,max_mean_se Minimum
#' @param prior_sample_size Prior sample size of strength.
#' @param min_prior_sample_size Minimum prior sample size or strength.
#' @param max_prior_sample_size Maximum prior sample size or strength.
#' @param slope_prior_sample_size Rate of change in prior
#' sample size or strength.
#' and maximum prior mean for sensitivity (se) and specificity (sp).
#' @importFrom magrittr %>%
#' @export
.get_prior_parameters <- function(thresholds,
                                  constant = TRUE,
                                  n_thresholds = NULL,
                                  n_models_or_tests = NULL,
                                  prior_p = NULL,
                                  prior_se = NULL,
                                  prior_sp = NULL,
                                  shift = 0.45,
                                  slope = 0.025,
                                  prior_sample_size = NULL,
                                  min_prior_sample_size = 20,
                                  max_prior_sample_size = 100,
                                  slope_prior_sample_size = 300,
                                  min_mean_se = 0.1,
                                  max_mean_se = 0.9,
                                  min_mean_sp = 0.1,
                                  max_mean_sp = 0.9) {
  if (isTRUE(constant)) {
    .priors <- .get_constant_prior_parameters(
      prior_p = prior_p,
      prior_se = prior_se,
      prior_sp = prior_sp,
      n_thresholds = length(thresholds),
      n_models_or_tests = n_models_or_tests
    )
  } else {
    .priors <- .get_nonconstant_prior_parameters(
      thresholds = thresholds,
      n_models_or_tests = n_models_or_tests,
      shift = shift,
      slope = slope,
      prior_sample_size = prior_sample_size,
      min_prior_sample_size = min_prior_sample_size,
      max_prior_sample_size = max_prior_sample_size,
      slope_prior_sample_size = slope_prior_sample_size,
      min_mean_se = min_mean_se,
      max_mean_se = max_mean_se,
      min_mean_sp = min_mean_sp,
      max_mean_sp = max_mean_sp
    )
  }

  return(.priors)
}

#' @title Get constant priors for Bayesian DCA
#'
#' @param n_thresholds Number of thresholds (int.).
#' @param n_models_or_tests Number of models or tests (int.).
#' @param prior_p,prior_se,prior_sp Non-negative shape values for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively.
#' Default is uniform prior for all parameters - Beta(1, 1).
#' A single vector of the form `c(a, b)` can be provided for each.
#' @importFrom magrittr %>%
#' @keywords internal
.get_constant_prior_parameters <- function(n_thresholds,
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


#' @title Get threshold- and model-specific priors for Bayesian DCA
#'
#' @param n_thresholds Number of thresholds (int.).
#' @param n_models_or_tests Number of models or tests (int.).
#' @param prior_p,prior_se,prior_sp Non-negative shape values for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively.
#' Default is uniform prior for all parameters - Beta(1, 1).
#' A single vector of the form `c(a, b)` can be provided for each.
#' @param min_mean_se,min_mean_sp,max_mean_se,max_mean_se Minimum
#' and maximum prior mean for sensitivity (se) and specificity (sp).
#' @param prior_sample_size Prior sample size or strength.
#' @param min_prior_sample_size Minimum prior sample size or strength.
#' @param max_prior_sample_size Maximum prior sample size or strength.
#' @param slope_prior_sample_size Rate of change in prior
#' sample size or strength.
#' @importFrom magrittr %>%
#' @keywords internal
.get_nonconstant_prior_parameters <- function(thresholds, # nolint
                                              n_models_or_tests,
                                              shift = 0.45,
                                              slope = 0.025,
                                              prior_sample_size = NULL,
                                              min_prior_sample_size = 20,
                                              max_prior_sample_size = 100,
                                              slope_prior_sample_size = 300,
                                              prior_p = NULL,
                                              min_mean_se = 0.1,
                                              max_mean_se = 0.9,
                                              min_mean_sp = 0.1,
                                              max_mean_sp = 0.9) {
  if (is.null(prior_p)) prior_p <- c(1, 1)

  stopifnot(
    "`shift` must be either a single number of a vector of size `n_models_or_tests`" = length(shift) == 1 | length(shift) == n_models_or_tests # nolint
  )
  stopifnot(
    "`slope` must be either a single number of a vector of size `n_models_or_tests`" = length(slope) == 1 | length(slope) == n_models_or_tests # nolint
  )
  stopifnot(
    "if given, `prior_sample_size` must be either a single number of a vector of size `n_models_or_tests`" = length(prior_sample_size) %in% c(0, 1, n_models_or_tests) # nolint
  )

  # if not model-specific, then use the first shift/slope
  if (length(shift) == 1L) {
    shift <- rep(shift, n_models_or_tests)
  }
  if (length(slope) == 1L) {
    slope <- rep(slope, n_models_or_tests)
  }
  if (is.null(prior_sample_size)) {
    prior_sample_size <- get_prior_sample_size(
      thresholds = thresholds,
      min_prior_sample_size = min_prior_sample_size,
      max_prior_sample_size = max_prior_sample_size,
      slope_prior_sample_size = slope_prior_sample_size
    )
  }
  if (length(prior_sample_size) == 1L) {
    prior_sample_size <- rep(prior_sample_size, n_models_or_tests)
  }
  n_thresholds <- length(thresholds)
  .priors <- list(
    p1 = prior_p[1],
    p2 = prior_p[2],
    Se1 = matrix(nrow = n_thresholds, ncol = n_models_or_tests),
    Se2 = matrix(nrow = n_thresholds, ncol = n_models_or_tests),
    Sp1 = matrix(nrow = n_thresholds, ncol = n_models_or_tests),
    Sp2 = matrix(nrow = n_thresholds, ncol = n_models_or_tests)
  )
  for (m in 1:n_models_or_tests) {
    se_mu <- get_prior_se_mu(
      thresholds = thresholds,
      shift = shift[m],
      slope = slope[m],
      .min = min_mean_se,
      .max = max_mean_se
    )
    sp_mu <- get_prior_sp_mu(
      thresholds = thresholds,
      shift = shift[m],
      slope = slope[m],
      .min = min_mean_sp,
      .max = max_mean_sp
    )
    smpl_size <- prior_sample_size[m]
    .priors[["Se1"]][, m] <- se_mu * smpl_size
    .priors[["Se2"]][, m] <- (1 - se_mu) * smpl_size
    .priors[["Sp1"]][, m] <- sp_mu * smpl_size
    .priors[["Sp2"]][, m] <- (1 - sp_mu) * smpl_size
  }

  return(.priors)
}
