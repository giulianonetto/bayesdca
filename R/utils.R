#' @title Compute Expected Value of Perfect Information (EVPI)
#'
#' @param obj BayesDCAList object
#' @param models_or_tests Character vector with models or tests to compare. If null, compares either first two in `obj$model_or_test_names` or the first one against Treat all/none (if only one available).
#' @importFrom magrittr %>%
evpi <- function(obj, models_or_tests = NULL) {
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

  .evpi <- vector("numeric", length(obj$thresholds))
  for (i in seq_along(obj$thresholds)) {
    # get posterior NB for each strategy in i-th threshold
    nb1 <- obj$draws$net_benefit[[models_or_tests[1]]][, i]
    if (length(models_or_tests) == 2) {
      nb2 <- obj$draws$net_benefit[[models_or_tests[2]]][, i]
    } else {
      nb2 <- obj$draws$treat_all[, i]
    }
    # EVPI equation in https://arxiv.org/abs/2208.03343 (p. 12)
    ENB_perfect <- mean(pmax(0, nb1, nb2))
    ENB_current <- max(0, mean(nb1), mean(nb2))
    .evpi[i] <- ENB_perfect - ENB_current
  }

  return(.evpi)
}


#' @title Plot Expected Value of Perfect Information (EVPI)
#'
#' @param obj BayesDCAList object
#' @param models_or_tests Character vector with models or tests to compare. If null, compares either first two in `obj$model_or_test_names` or the first one against Treat all/none (if only one available).
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
plot_evpi <- function(obj, models_or_tests = NULL, labels = NULL) {
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

  # get subtitles
  if (length(models_or_tests) == 1) {
    .subtitle <- paste0("EVPI: ",
                        plot_labels,
                        ' vs. Treat all or none')
  } else {
    .subtitle <- paste0("EVPI: ",
                        plot_labels[1],
                        ' vs. ',
                        plot_labels[2])
  }

  data.frame(
    .threhsolds = obj$thresholds,
    .evpi = evpi(obj, models_or_tests = models_or_tests)
  ) %>%
    ggplot2::ggplot(ggplot2::aes(.threhsolds, .evpi)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(10)
    ) +
    ggplot2::labs(
      x = "Decision threshold", y = NULL,
      subtitle = .subtitle
    )
}


#' @title Get colors and labels for BayesDCA plots
#'
#' @param obj BayesDCAList object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
get_colors_and_labels <- function(obj, models_or_tests = NULL, colors = NULL, labels = NULL) {
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
  color_values <- c(
    "Treat all" = "black", "Treat none" = "gray40"
  )
  n_colors <- length(model_or_test_names)
  if (n_colors < 9) {
    palette <- RColorBrewer:::brewer.pal(max(c(n_colors, 3)), 'Dark2')
  } else {
    palette <- grDevices::colorRampPalette(
      RColorBrewer:::brewer.pal(n_colors, 'Set2')
    )(n_colors)
  }
  # set actual color values to use in scale_color_manual
  for (i in seq_len(n_colors)) {
    model_or_test <- model_or_test_names[i]
    if (!is.null(colors) & model_or_test %in% names(colors)) {
      color_values[[model_or_test]] <- colors[[model_or_test]]
    } else {
      color_values[[model_or_test]] <- palette[i]
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
get_survival_time_exposed <- function(.prediction_time, .cutpoints) {
  time_exposed <- numeric(length = length(cutpoints))
  for (i in 1:(length(cutpoints)-1) ) {
    if (prediction_time >= cutpoints[i+1]) {
      # if prediction time > upper bound, use interval length
      time_exposed[i] <- cutpoints[i+1] - cutpoints[i]
    } else if (prediction_time >= cutpoints[i]) {
      # if prediction time > lower bound, use distance between pred time and lower bound
      time_exposed[i] <- prediction_time - cutpoints[i]
    } else {
      # otherwise, prediction time is before interval, exposure time is zero
      time_exposed[i] <- 0
    }
  }
  return(time_exposed)
}

#' @title Get posterior parameters for survival model
#'
#' @param .prediction_data Contains prognostic model predictions
#' @param .surv_data Contains observed time to event (`.time` column) and event indicator (`.status` column)
#' @param .cutpoints Survival cutpoints.
#' @param .models_or_tests Names for models or tests.
#' @param .thresholds DCA thresholds.
#' @param .prior_scaling_factor Prior alpha = (0.69/median_surv)*.prior_scaling_factor, Prior beta = .prior_scaling_factor.
#' @importFrom magrittr %>%
#' @importFrom survival Surv
get_survival_posterior_parameters <- function(
    .prediction_data,
    .surv_data,
    .cutpoints,
    .models_or_tests,
    .thresholds,
    .prior_scaling_factor = 0.1
) {

  initialize_pars <- function() {
    lapply(
      seq_along(.models_or_tests),
      function(...) {
        matrix(nrow = length(.cutpoints),
               ncol = length(.thresholds))
      }
    )
  }

  all_posterior_alphas <- initialize_pars()
  all_posterior_betas <- initialize_pars()

  for (i in seq_along(.models_or_tests)) {

    .model <- .models_or_tests[i]

    for (j in seq_along(.thresholds)) {
      .thr <- .thresholds[j]
      .predictions <- .prediction_data[[.model]]
      .positive_prediction <- .predictions >= .thr
      .d <- .surv_data[.positive_prediction, ]
      .d$patient_id <- 1:nrow(.d)
      .median_survival <- survival:::median.Surv(Surv(.d$.time, .d$.status))$quantile
      .prior_alpha <- (0.69/.median_survival) * .prior_scaling_factor
      .prior_beta <- .prior_scaling_factor
      .d_split <- survival::survSplit(
        Surv(.time, .status) ~ 1,
        data = .d,
        cut = .cutpoints,
        subset = 1:nrow(.d),
        id = "patient_id",
        start = "tstart",
        end = "tstop"
      ) %>%
        # TODO: might want to double check the calculation below
        dplyr::group_by(
          interval_id = paste0(
            "interval_", as.numeric(factor(tstart))
          )
        ) %>%
        dplyr::summarise(
          total_events = sum(.status),
          total_exposure_time = sum(tstop - tstart)
        )

      all_posterior_alphas[[i]][ , j] <- .d_split$total_events + .prior_alpha
      all_posterior_betas[[i]][ , j] <- .d_split$total_events + .prior_beta

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
get_positivity_posterior_parameters <- function(
    .prediction_data,
    .thresholds,
    .prior_shape1 = 0.5,
    .prior_shape2 = 0.5
) {

  N <- nrow(.prediction_data)
  n_models <- ncol(.prediction_data)
  .models_or_tests <- colnames(.prediction_data)
  n_thresholds <- length(.thresholds)

  all_posterior_shape1 <- matrix(nrow = n_models,
                                 ncol = n_thresholds)
  all_posterior_shape2 <- matrix(nrow = n_models,
                                 ncol = n_thresholds)

  for (i in seq_along(.models_or_tests)) {
    .model <- .models_or_tests[i]
    for (j in seq_along(.thresholds)) {
      .thr <- .thresholds[j]
      .predictions <- .prediction_data[[.model]]
      .positive_prediction <- .predictions >= .thr
      total_positives <- sum(.positive_prediction)
      all_posterior_shape1[i, j] <- total_positives + .prior_shape1
      all_posterior_shape2[i, j] <- N - total_positives + .prior_shape2
    }
  }

  .posterior_pars <- list(
    .shape1 = all_posterior_shape1,
    .shape2 = all_posterior_shape2
  )
  return(.posterior_pars)
}
