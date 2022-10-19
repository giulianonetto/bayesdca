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
