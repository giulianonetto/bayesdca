#' @title Plot BayesDCA
#'
#' @param obj BayesDCA object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @param models_or_tests Character vector with
#' models or tests to compare. If null, compares
#' either first two in `obj$model_or_tests`
#' or the first one against Treat all/none
#' (if only one available).
#' @importFrom magrittr %>%
#' @export
plot.BayesDCA <- function(obj,
                          models_or_tests = NULL,
                          colors = NULL,
                          labels = NULL,
                          raw_values = NULL,
                          raw_values_label = "Biomarker threshold",
                          linewidth = 1.5, ...) {
  models_or_tests <- validate_models_or_tests(
    obj = obj, models_or_tests = models_or_tests
  )

  net_benefit_data <- obj$summary$net_benefit %>%
    dplyr::filter(model_or_test_name %in% models_or_tests)

  colors_and_labels <- get_colors_and_labels(
    obj = obj,
    colors = colors,
    labels = labels,
    models_or_tests = models_or_tests
  )

  .ymin <- ifelse(
    max(obj$summary$treat_all$estimate) > 0.02,
    -0.02,
    -max(obj$summary$treat_all$estimate)
  )

  .p <- ggplot2::ggplot() +
    # set x axis
    ggplot2::aes(x = threshold) + # nolint
    # add color/fill/label scheme
    colors_and_labels +
    # add treat all curve
    ggplot2::geom_ribbon(
      data = obj$summary$treat_all,
      ggplot2::aes(
        ymax = `97.5%`, ymin = `2.5%`, # nolint
        fill = "Treat all"
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = obj$summary$treat_all, linewidth = linewidth,
      ggplot2::aes(y = estimate, color = "Treat all", group = 1) # nolint
    ) +
    # add net benefit curves
    ggplot2::geom_ribbon(
      data = net_benefit_data,
      ggplot2::aes(
        ymin = `2.5%`, ymax = `97.5%`,
        fill = model_or_test_name # nolint
      ), # nolint
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = net_benefit_data,
      linewidth = linewidth,
      ggplot2::aes(
        y = estimate,
        color = model_or_test_name,
        group = model_or_test_name
      )
    ) +
    # add treat none curve
    ggplot2::geom_hline(
      ggplot2::aes(color = "Treat none", yintercept = 0),
      linetype = 2, linewidth = linewidth,
    ) +
    # make it pretty
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::coord_cartesian(ylim = c(.ymin, NA)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(
      x = "Decision threshold", y = "Net Benefit",
      color = NULL
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        keywidth = ggplot2::unit(1, "cm"),
        override.aes = list(linewidth = 2)
      )
    )

  if (!is.null(raw_values)) {
    stopifnot(is.data.frame(raw_values))

    raw_values_p <- ggplot2::ggplot(
      raw_values,
      ggplot2::aes(x = thresholds * 100) # nolint
    ) +
      ggplot2::scale_x_continuous(
        labels = raw_values$values,
      ) +
      ggplot2::theme_bw(base_size = 14) +
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

#' @title Plot BayesDCASurv
#'
#' @param obj BayesDCASurv object
#' @export
plot.BayesDCASurv <- function(obj,
                              models_or_tests = NULL,
                              colors = NULL,
                              labels = NULL,
                              raw_values = NULL,
                              raw_values_label = "Biomarker threshold",
                              linewidth = 1.5) {
  models_or_tests <- validate_models_or_tests(
    obj = obj, models_or_tests = models_or_tests
  )

  net_benefit_data <- obj$summary$net_benefit %>%
    dplyr::filter(model_or_test_name %in% models_or_tests)

  colors_and_labels <- get_colors_and_labels(
    obj = obj,
    colors = colors,
    labels = labels,
    models_or_tests = models_or_tests
  )

  .ymin <- ifelse(
    max(obj$summary$treat_all$estimate) > 0.02,
    -0.02,
    -max(obj$summary$treat_all$estimate)
  )

  .p <- ggplot2::ggplot() +
    # set x axis
    ggplot2::aes(x = threshold) + # nolint
    # add color/fill/label scheme
    colors_and_labels +
    # add treat all curve
    ggplot2::geom_ribbon(
      data = obj$summary$treat_all,
      ggplot2::aes(
        ymax = `97.5%`, ymin = `2.5%`, # nolint
        fill = "Treat all"
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = obj$summary$treat_all, linewidth = linewidth,
      ggplot2::aes(y = estimate, color = "Treat all", group = 1) # nolint
    ) +
    # add net benefit curves
    ggplot2::geom_ribbon(
      data = net_benefit_data,
      ggplot2::aes(
        ymin = `2.5%`, ymax = `97.5%`,
        fill = model_or_test_name # nolint
      ), # nolint
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = net_benefit_data,
      linewidth = linewidth,
      ggplot2::aes(
        y = estimate,
        color = model_or_test_name,
        group = model_or_test_name
      )
    ) +
    # add treat none curve
    ggplot2::geom_hline(
      ggplot2::aes(color = "Treat none", yintercept = 0),
      linetype = 2, linewidth = linewidth,
    ) +
    # make it pretty
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::coord_cartesian(ylim = c(.ymin, NA)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(
      x = "Decision threshold", y = "Net Benefit",
      color = NULL
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        keywidth = ggplot2::unit(1, "cm"),
        override.aes = list(linewidth = 2)
      )
    )

  if (!is.null(raw_values)) {
    stopifnot(is.data.frame(raw_values))

    raw_values_p <- ggplot2::ggplot(
      raw_values,
      ggplot2::aes(x = thresholds * 100) # nolint
    ) +
      ggplot2::scale_x_continuous(
        labels = raw_values$values,
      ) +
      ggplot2::theme_bw(base_size = 14) +
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

#' @title Plot BayesDCA comparison
#'
#' @param obj BayesDCA object
#' @param models_or_tests Character vector with models or tests to compare. If null, compares either first two in `obj$model_or_tests` or the first one against Treat all/none (if only one available).
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
compare_dca <- function(obj,
                        models_or_tests = NULL,
                        colors = NULL,
                        labels = NULL,
                        plot_list = FALSE,
                        .evpi = FALSE,
                        type = c("best", "useful", "pairwise"),
                        ...) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two models_or_tests to plot pairwise comparison" = length(models_or_tests) == 2
    )
  }

  if (is.null(models_or_tests)) {
    models_or_tests <- as.vector(na.omit(obj$model_or_tests))
  } else {
    stopifnot(
      "Provided `models_or_tests` are not available" = all(
        models_or_tests %in% obj$model_or_tests
      )
    )
  }



  p1 <- plot(obj,
    colors = colors, labels = labels,
    models_or_tests = models_or_tests
  )
  p2 <- plot_delta(
    obj = obj,
    models_or_tests = models_or_tests,
    colors = colors,
    labels = labels,
    type = type
  ) +
    ggplot2::guides(color = "none", fill = "none")
  p3 <- plot_superiority_prob(
    obj = obj,
    models_or_tests = models_or_tests,
    colors = colors,
    labels = labels,
    type = type
  ) +
    ggplot2::guides(color = "none")

  if (isTRUE(plot_list)) {
    .plot_list <- list(
      dca = p1,
      delta = p2,
      prob_better = p3
    )
    if (isTRUE(.evpi)) {
      .plot_list[["evpi"]] <- plot_evpi(
        obj = obj,
        models_or_tests = models_or_tests,
        colors = colors,
        labels = labels,
        type = type
      )
    }

    return(.plot_list)
  } else {
    if (isTRUE(.evpi)) {
      p1 <- p1 + ggplot2::labs(subtitle = "DCA")
      p4 <- plot_evpi(
        obj = obj,
        models_or_tests = models_or_tests,
        colors = colors,
        labels = labels,
        type = type
      ) +
        ggplot2::guides(color = "none")
      .p <- (p1 | p3) / (p2 | p4) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(
          legend.position = "bottom",
          legend.text = ggplot2::element_text(size = 14)
        )
    } else {
      .p <- p1 / (p2 | p3) +
        patchwork::plot_layout(guides = "collect")
    }
    return(.p)
  }
}

#' @title Plot BayesDCA delta
#'
#' @param obj BayesDCA object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot_delta(fit)
#' @return A ggplot object.
plot_delta <- function(obj,
                       models_or_tests = NULL,
                       type = c("best", "useful", "pairwise"),
                       colors = NULL,
                       labels = NULL) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two models_or_tests to plot pairwise comparison" = length(models_or_tests) == 2 # nolint
    )
  }

  models_or_tests <- validate_models_or_tests(
    obj = obj, models_or_tests = models_or_tests
  )

  if (is.null(labels)) {
    labels <- setNames(models_or_tests, models_or_tests)
  } else {
    stopifnot(
      "Names of labels must match models_or_tests" = all(sort(names(labels)) == sort(models_or_tests)) # nolint
    )
  }

  if (inherits(obj, "BayesDCA")) {
    df <- get_delta_plot_data_binary(
      obj = obj,
      models_or_tests = models_or_tests,
      type = type,
      labels = labels
    )
  } else if (inherits(obj, "BayesDCASurv")) {
    df <- get_delta_plot_data_surv(
      obj = obj,
      models_or_tests = models_or_tests,
      type = type,
      labels = labels
    )
  } else {
    msg <- paste0(
      "FATAL - unknown object: ", class(obj),
      "\nIt should be either 'BayesDCA' or 'BayesDCASurv'."
    )
    stop(msg)
  }

  if (type == "pairwise") {
    .subtitle <- paste0(
      labels[models_or_tests[1]],
      " v.s. ",
      labels[models_or_tests[2]]
    )
    initial_plot <- df %>%
      ggplot2::ggplot() +
      ggplot2::aes(
        x = threshold, y = estimate, # nolint
        ymin = `2.5%`, ymax = `97.5%` # nolint
      ) + # nolint
      ggplot2::geom_ribbon(alpha = 0.3) +
      ggplot2::geom_line(lwd = 0.9)
  } else {
    .subtitle <- paste0(
      "Difference against ",
      ifelse(
        type == "useful",
        "treat all or none",
        "best competitor"
      )
    )
    .colors_and_labels <- get_colors_and_labels(
      obj = obj,
      colors = colors,
      labels = labels,
      models_or_tests = models_or_tests,
      all_or_none = FALSE
    )
    initial_plot <- df %>%
      ggplot2::ggplot() +
      ggplot2::aes(
        x = threshold, y = estimate, # nolint
        ymin = `2.5%`, ymax = `97.5%` # nolint
      ) + # nolint
      ggplot2::geom_ribbon(alpha = 0.3, ggplot2::aes(fill = model_or_test)) + # nolint
      ggplot2::geom_line(ggplot2::aes(color = model_or_test), lwd = 0.9) +
      .colors_and_labels
  }

  .ymax <- max(
    obj$summary$treat_all$estimate
  ) * 1.02 # put into context of max treat all
  .ymin <- -1 * .ymax

  .plot <- initial_plot +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 2,
      color = "gray40", lwd = 0.8
    ) +
    ggplot2::labs(
      x = "Decision threshold",
      y = expression(Delta[NB]),
      subtitle = .subtitle,
      color = NULL, fill = NULL
    ) +
    ggplot2::coord_cartesian(
      ylim = c(.ymin, .ymax)
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        keywidth = ggplot2::unit(1, "cm"),
        override.aes = list(linewidth = 2)
      )
    )
  return(.plot)
}

#' @title Plot P(useful) from Wynants 2018 (doi: 10.1002/sim.7653)
#' @param obj BayesDCA or BayesDCASurv object
#' @param min_diff Minimal difference for superiority. Defaults to zero. Used only for `type = "pairwise"`
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(PredModelData)
#' fit <- dca(PredModelData, cores = 4)
#' plot_superiority_prob(fit)

#' @return A ggplot object.
plot_superiority_prob <- function(obj, models_or_tests = NULL, type = c("best", "useful", "pairwise"), min_diff = 0, colors = NULL, labels = NULL) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two models_or_tests to plot pairwise comparison" = length(models_or_tests) == 2 # nolint
    )
  }

  models_or_tests <- validate_models_or_tests(
    obj = obj, models_or_tests = models_or_tests
  )

  if (is.null(labels)) {
    labels <- setNames(models_or_tests, models_or_tests)
  } else {
    stopifnot(
      "Names of labels must match models_or_tests" = all(names(labels) == models_or_tests)
    )
  }

  if (inherits(obj, "BayesDCA")) {
    df <- get_superiority_prob_plot_data_binary(
      obj = obj,
      min_diff = min_diff,
      models_or_tests = models_or_tests,
      type = type,
      labels = labels
    )
  } else if (inherits(obj, "BayesDCASurv")) {
    df <- get_superiority_prob_plot_data_surv(
      obj = obj,
      models_or_tests = models_or_tests,
      type = type,
      labels = labels
    )
  } else {
    msg <- paste0(
      "FATAL - unknown object: ", class(obj),
      "\nIt should be either 'BayesDCA' or 'BayesDCASurv'."
    )
    stop(msg)
  }



  if (type == "pairwise") {
    .subtitle <- paste0("P(", labels[models_or_tests[1]], " better than ", labels[models_or_tests[2]], ")") # nolint
    initial_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = threshold, y = prob)) + # nolint
      ggplot2::geom_line(lwd = 0.9)
  } else {
    .subtitle <- paste0("P(", type, ")")
    .colors_and_labels <- get_colors_and_labels(
      obj = obj,
      colors = colors,
      labels = labels,
      models_or_tests = models_or_tests,
      all_or_none = FALSE
    )
    initial_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = threshold, y = prob)) + # nolint
      ggplot2::geom_line(ggplot2::aes(color = model_or_test), lwd = 0.9) + # nolint
      .colors_and_labels
  }


  .plot <- initial_plot +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(1),
      breaks = seq(0, 1, 0.1),
      limits = c(0, 1.01)
    ) +
    ggplot2::labs(
      x = "Decision threshold",
      y = NULL,
      color = NULL,
      subtitle = .subtitle
    )

  return(.plot)
}

#' @title Plot Expected Value of Perfect Information (EVPI)
#'
#' @param obj BayesDCA object
#' @param models_or_tests Character vector with models or tests
#' to compare. If null, compares either first two in
#' `obj$model_or_tests` or the first one against
#' Treat all/none (if only one available).
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
plot_evpi <- function(obj, models_or_tests = NULL, type = c("best", "useful", "pairwise"), colors = NULL, labels = NULL) { # nolint
  type <- match.arg(type)

  if (type != "best") {
    msg <- paste0(
      "CAREFUL: EVPI might only make sense for type='best'. Use type='",
      type,
      "' at your own risk."
    )
    message(cli::col_br_cyan(msg))
  }
  if (type == "pairwise") {
    stopifnot(
      "Must specify two models_or_tests to plot pairwise comparison" = length(models_or_tests) == 2 # nolint
    )
  }

  if (is.null(obj$draws)) {
    msg <- "Retrieving posterior draws."
    message(msg)
    if (inherits(obj, "BayesDCA")) {
      obj$draws <- .extract_dca_draws(
        fit = obj,
        model_or_tests = models_or_tests
      )
    } else {
      obj$draws <- .extract_dca_surv_draws(
        fit = obj,
        model_or_tests = models_or_tests
      )
    }
  }

  if (is.null(models_or_tests)) {
    models_or_tests <- as.vector(na.omit(obj$model_or_tests))
  } else {
    stopifnot(
      "Provided `models_or_tests` are not available" = all(
        models_or_tests %in% obj$model_or_tests
      )
    )
  }

  if (is.null(labels)) {
    labels <- setNames(models_or_tests, models_or_tests)
  } else {
    stopifnot(
      "Names of labels must match models_or_tests" = all(names(labels) == models_or_tests)
    )
  }

  if (type == "pairwise") {
    nb1 <- obj$draws$net_benefit[[models_or_tests[1]]]
    nb2 <- obj$draws$net_benefit[[models_or_tests[2]]]
    df <- tibble::tibble(
      .evpi = evpi(thresholds = obj$thresholds, nb1, nb2),
      threshold = obj$thresholds
    )
    .subtitle <- paste0("EVPI: ", labels[models_or_tests[1]], " v.s. ", labels[models_or_tests[2]]) # nolint
    initial_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = threshold, y = .evpi)) + # nolint
      ggplot2::geom_line(lwd = 0.9)
  } else if (type == "useful") {
    df <- lapply(
      seq_along(models_or_tests),
      function(i) {
        .m <- models_or_tests[i]

        # type = 'useful', considers this model against treat all/none
        args <- list(
          thresholds = obj$thresholds,
          obj$draws$net_benefit[[.m]],
          obj$draws$treat_all
        )

        .evpi <- do.call(evpi, args)


        tibble::tibble(
          .evpi = .evpi,
          threshold = obj$thresholds,
          model_or_test = .m,
          label = labels[models_or_tests[i]]
        )
      }
    ) %>%
      dplyr::bind_rows()

    .colors_and_labels <- get_colors_and_labels(
      obj = obj,
      colors = colors,
      labels = labels,
      models_or_tests = models_or_tests,
      all_or_none = FALSE
    )
    .subtitle <- paste0("EVPI against treat all or none")
    initial_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = threshold, y = .evpi)) + # nolint
      ggplot2::geom_line(ggplot2::aes(color = model_or_test), lwd = 0.9) + # nolint
      .colors_and_labels
  } else {
    # actual EVPI
    args <- obj$draws$net_benefit
    args[[length(obj$draws$net_benefit) + 1]] <- obj$draws$treat_all
    args[["thresholds"]] <- obj$thresholds
    df <- tibble::tibble(
      .evpi = do.call(evpi, args),
      threshold = obj$thresholds,
      model_or_test = NA_character_,
      label = NA_character_
    )
    .subtitle <- "EVPI"
    initial_plot <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = threshold, y = .evpi)) + # nolint
      ggplot2::geom_line(lwd = 0.9)
  }



  .plot <- initial_plot +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(10)
    ) +
    ggplot2::labs(
      x = "Decision threshold",
      y = NULL,
      color = NULL,
      subtitle = .subtitle
    )

  return(.plot)
}

#' Plot classification performance from Bayesian DCA of binary outcome
#'
#' May plot either sensitivity or specificity.
#' @param obj BayesDCA object
#' @param colors Named vector with color for each model or test. If provided
#' for a subset of models or tests, only that subset will be plotted.
#' @param labels Named vector with label for each model or test.
#' @param models_or_tests Character vector with models or tests
#' to compare. If null, compares either first two in
#' `obj$model_or_tests` or the first one against
#' Treat all/none (if only one available).
#' @importFrom magrittr %>%
#' @export
plot_classification <- function(obj,
                                type = c("sensitivity", "specificity"),
                                models_or_tests = NULL,
                                colors = NULL,
                                labels = NULL) {
  type <- match.arg(type)

  if (!is.null(models_or_tests)) {
    stopifnot(
      "Provided `models_or_tests` are not available" = all(
        models_or_tests %in% obj$model_or_tests
      )
    )

    plot_data <- obj$summary[[type]] %>%
      dplyr::filter(model_or_test_name %in% models_or_tests) # nolint
  } else {
    models_or_tests <- obj$model_or_tests
    plot_data <- obj$summary[[type]]
  }

  colors_and_labels <- get_colors_and_labels(
    obj = obj,
    colors = colors,
    labels = labels,
    models_or_tests = models_or_tests,
    all_or_none = FALSE
  )
  .p <- plot_data %>%
    dplyr::filter(
      model_or_test_name %in% models_or_tests # nolint
    ) %>%
    ggplot2::ggplot() +
    # set x axis
    ggplot2::aes(x = threshold) + # nolint
    # add color/fill/label scheme
    colors_and_labels +
    # add net benefit curves
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = `2.5%`, ymax = `97.5%`, # nolint
        fill = model_or_test_name
      ),
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = estimate, # nolint
        color = model_or_test_name,
        group = model_or_test_name
      )
    ) +
    # make it pretty
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::labs(
      x = "Decision threshold",
      y = stringr::str_to_title(type),
      color = NULL
    ) +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(
        keywidth = ggplot2::unit(0.5, "cm")
      )
    )

  return(.p)
}
