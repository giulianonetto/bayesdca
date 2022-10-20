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

get_survival_cutpoints <- function(death_times) {
  if (length(death_times) >= 200) {
    .probs <- c(.1, .25, .5, .75, .9, .975)
    
    top_one_percent <- quantile(death_times, 0.99)
    if (sum(death_times >= top_one_percent) > 50) {
      .probs <- c(.probs, .99)
    } 
  } else if (length(death_times) >= 50) {
    .probs <- c(.1, .5, .9)
  } else if (length(death_times) >= 30) {
    .probs <- c(1/3, 2/3)
  } else {
    .probs <- c(1/3)
  }
  
  .cuts <- c(
    0,
    unname(quantile(death_times, probs = .probs))
  )
  
  return(.cuts)
}


get_survival_time_exposed <- function(.prediction_time, .cutpoints) {
  # reference: https://rpubs.com/kaz_yos/surv_stan_piecewise1
  # TODO: <= or <
  if (!(0 %in% .cutpoints)) {
    .cutpoints <- c(0, .cutpoints)
  }
  interval_exposed <- outer(.cutpoints, .prediction_time, `<=`)
  
  ## t - cutpoint. Multiply by interval exposed to avoid negative times.
  time_exposed <-  -outer(.cutpoints, .prediction_time, `-`) * interval_exposed
  if (any(is.na(time_exposed))) {
    cat("There are NAs in time exposed. Replacing with zero.\n")
    time_exposed[is.na(time_exposed)] <- 0.0
  }
  
  
  ## Last interval is of width Inf
  interval_widths <- c(diff(.cutpoints), Inf)
  
  ## For each interval, time exposed cannot exceed interval width.
  ## matrix with length(.cutpoints) rows and length(.t) columns
  time_exposed_correct  <- sweep(x = time_exposed,
                                 MARGIN = 1,
                                 STATS = interval_widths,
                                 FUN = pmin)
  
  return(t(time_exposed_correct))
}

get_survival_posterior_parameters <- function(
    .data, 
    .cutpoints,
    .models_or_tests,
    .thresholds,
    .prior_alpha, 
    .prior_beta
) {
  
  all_posteriorpars <- vector('list', length(.models_or_tests))
  for (i in seq_along(.models_or_tests)) {
    for (j in seq_along(.thresholds)) {
      .model <- .models_or_tests[i]
      .thr <- .thresholds[j]
      .predictions <- .data[[.model]]
      .d <- .data[.predictions >= .thr, ]
      .d$patient_id <- 1:nrow(.d)
      .median_surv <- median(.d[["outcomes"]])$quantile
      .d_split <- survSplit(
        outcomes ~ 1,
        data = .d,
        cut = .cutpoints,
        subset = 1:nrow(.d),
        id = "patient_id",
        start = "tstart",
        end = "tstop"
      )
      # TODO: group_by and so on
    }
  }
  lapply(.thresholds, \(i) {
    .d <- .data %>% 
      dplyr::filter(cancerpredmarker >= i)
    q50_surv <- median(.d$outcomes)$quantile
    survSplit(
      s ~ 1,
      data = .d,
      cut = .cutpoints,
      subset = 1:nrow(.d),
      id = "patientid",
      start = "tstart",
      end = "tstop"
    ) %>% 
      group_by(
        ij = paste0("[", tstart, ", ", tstop, ")"),
        j = paste0(
          "interval_", as.numeric(factor(tstart))
        )
      ) %>% 
      summarise(
        dij = sum(cancer),
        tij = sum(tstop - tstart),
        .groups = 'drop'
      ) %>% 
      ungroup() %>% 
      group_by(j) %>% 
      summarise(
        total_dij = sum(dij),
        total_tij = sum(tij)
      ) %>% 
      mutate(.thr = i)
    
    posterior_parameters <- data.frame(
      alpha = .prior_alpha + .totals$total_dij,
      beta = .prior_beta + .totals$total_tij
    )
    
    return(posterior_parameters)
  })
  
  
}
