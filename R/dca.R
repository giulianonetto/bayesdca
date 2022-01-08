#' Fit Bayesian Decision Curve Analysis using Stan
#'
#' @param N Integer specifying total number of samples (i.e., participants).
#' @param d Diseased: integer specifying total number of diseased persons or events.
#' @param tp True Positives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model (positive count among cases).
#' @param tn True Negatives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model (negative count among non-cases).
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param model_type Character indicating model type. Either "ic"
#' ("independent conjugate", default) or "mo" ("multiple observations").
#' @param prior_p,prior_se,prior_sp Numeric vectors with shapes for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively. Default is c(1, 1) for all (uniform prior).
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#'
.dca_stan <- function(N, d, tp, tn, B = NULL,
                      thresholds = seq(0, 0.5, 0.01),
                      model_type = "ic",
                      prior_p = c(1, 1),
                      prior_se = c(1, 1),
                      prior_sp = c(1, 1),
                      refresh = 0, ...) {
  
  if (is.null(thresholds)) {
    thresholds <- 0
  }
  
  thresholds <- pmin(thresholds, 0.999)  # odds(1) = Inf
  
  n_thresholds <- length(thresholds)
  
  standata <- list(
    # input data
    N = N, d = d, tp = tp, tn = tn, B = B,
    # priors
    prior_p1 = prior_p[[1]], prior_p2 = prior_p[[2]],
    prior_Se1 = prior_se[[1]], prior_Se2 = prior_se[[2]],
    prior_Sp1 = prior_sp[[1]], prior_Sp2 = prior_sp[[2]],
    # info for generated quantities
    n_thresholds = n_thresholds,
    thresholds = array(thresholds)
  )
  
  if (model_type %in% c("correlated", "c")) {
    .model <- stanmodels$dca_correlated
  } else if (model_type %in% c("independent", "i")) {
    .model <- stanmodels$dca_independent
  } else if (model_type %in% c("independent conjugate", "ic")) {
    .model <- stanmodels$dca_independent_conjugate
  } else if (model_type %in% c("independent conjugate rates", "icr")) {
    .model <- stanmodels$dca_independent_conjugate_rates
  } else if (model_type %in% c("multiple obs", "mo")) {
    .model <- stanmodels$dca_independent_conjugate_multiple_obs
  } else {
    stop("Invalid `model_type`.")
  }
  
  stanfit <- rstan::sampling(.model,
                             data = standata,
                             refresh = refresh, ...)
  return(stanfit)
}

#' Bayesian Decision Curve Analysis for Binary Test
#'
#' @export
#' @param N Integer specifying total number of samples (i.e., number of participants).
#' @param d Diseased: integer specifying total number of diseased persons or outcome cases.
#' @param tp True Positives: integer specifying total number of diseased persons or cases correctly
#' identified as such by the binary test (positive tests among cases). `tp` must be at most `d`.
#' @param tn True Negatives: integer specifying total number of non-diseased persons or non-cases correctly
#' identified as such by the diagnostic test or prediction model (negative tests among non-cases). `tn` must be at most `N - d`.
#' @param thresholds Numeric vector with probability thresholds with which
#' the decision curve should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) object in
#' the output (default FALSE).
#' @param prior_p,prior_se,prior_sp Numeric vector with shapes for
#' Beta(alpha, beta) prior used for p, Se, and Sp. Default is c(1, 1) for all.
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (check its help
#' page for details).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class `DiagTestDCA`
#' @importFrom magrittr %>%
#' @examples
#'
#' fit <- dca_binary_test(N = 500, d = 83, tp = 77, tn = 378)
#' plot(fit)
dca_binary_test <- function(N, d, tp, tn,
                            thresholds = seq(0, 0.5, 0.01),
                            keep_fit = FALSE,
                            keep_draws = TRUE,
                            prior_p = c(1, 1),
                            prior_se = c(1, 1),
                            prior_sp = c(1, 1),
                            summary_probs = c(0.025, 0.975),
                            refresh = 0, ...) {
  
  thresholds <- pmin(thresholds, 0.99)
  
  fit <- .dca_stan(
    N = N, d = d, tp = tp, tn = tn, model_type = "ic",
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    thresholds = thresholds, refresh = refresh,
    ...
  )
  
  fit_summary <- rstan::summary(
    fit, probs = summary_probs
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("par_name") %>%
    tibble::as_tibble() %>%
    dplyr::select(-se_mean) %>%
    dplyr::rename(estimate := mean)
  
  model_parameters <- fit_summary %>%
    dplyr::filter(
      stringr::str_detect(par_name, "p|Sp|Se")
    )
  net_benefit <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "net_benefit")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  treat_all <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "treat_all")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  
  output_data <- list(
    model_parameters = model_parameters,
    net_benefit = net_benefit, treat_all = treat_all,
    N = N, d = d, tp = tp, tn = tn,
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    thresholds = thresholds
  )
  
  if (isTRUE(keep_fit)) {
    output_data[['fit']] <- fit
  }
  
  if(isTRUE(keep_draws)) {
    output_data[['draws']] <- rstan::extract(fit)
  }
  
  .output <- structure(output_data, class = "DiagTestDCA")
  return(.output)
}

#' Bayesian Decision Curve Analysis for Predictive Model
#'
#' @export
#' @param outcomes Integer vector with outcome events for each individual (0 or 1).
#' @param predictions Numeric vector with predicted probabilities.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default FALSE).
#' @param prior_p,prior_se,prior_sp Numeric vector with shapes for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively. Default is c(1, 1) for all (uniform prior).
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `DiagTestDCA`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' head(PredModelData)
#' fit <- dca_predictive_model(outcomes = PredModelData$outcomes,
#'                             predictions = PredModelData$predictions)
#' plot(fit)
dca_predictive_model <- function(outcomes,
                                 predictions,
                                 keep_draws = TRUE,
                                 keep_fit = FALSE,
                                 prior_p = c(1, 1),
                                 prior_se = c(1, 1),
                                 prior_sp = c(1, 1),
                                 summary_probs = c(0.025, 0.975),
                                 refresh = 0, ...) {
  
  df <- get_thr_data(outcomes = outcomes,
                     predictions = predictions)
  f <- function(i, n) lapply(i, function(k) rep(k, n))
  prior_se <- f(prior_se, nrow(df))
  prior_sp <- f(prior_sp, nrow(df))
  
  fit <- .dca_stan(
    N = df$N, d = df$d, tp = df$tp, tn = df$tn, B = nrow(df),
    thresholds = df$thresholds, model_type = "mo",
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    refresh = refresh,
    ...
  )
  
  fit_summary <- rstan::summary(
    fit, probs = summary_probs
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("par_name") %>%
    tibble::as_tibble() %>%
    dplyr::select(-se_mean) %>%
    dplyr::rename(estimate := mean)
  
  model_parameters <- fit_summary %>%
    dplyr::filter(
      stringr::str_detect(par_name, "p|Sp|Se")
    )
  net_benefit <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "net_benefit")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = df$thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  treat_all <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "treat_all")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = df$thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  
  output_data <- list(
    model_parameters = model_parameters,
    net_benefit = net_benefit, treat_all = treat_all,
    N = df$N, d = df$d, tp = df$tp, tn = df$tn,
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    thresholds = df$thresholds
  )
  
  if (isTRUE(keep_fit)) {
    output_data[['fit']] <- fit
  }
  
  if(isTRUE(keep_draws)) {
    output_data[['draws']] <- rstan::extract(fit)
  }
  
  .output <- structure(output_data, class = "PredModelDCA")
  return(.output)
}

#' Bayesian Decision Curve Analysis for a pair of Predictive Models
#'
#' @export
#' @param outcomes Integer vector with outcome events for each individual (0 or 1).
#' @param predictions1 Numeric vector with predicted probabilities for model 1.
#' @param predictions2 Numeric vector with predicted probabilities for model 2.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default FALSE).
#' @param prior_p,prior_se,prior_sp Numeric vector with shapes for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively. Default is c(1, 1) for all (uniform prior).
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `DiagTestDCA`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' head(PredModelData)
#' fit <- dca_predictive_model(outcomes = PredModelData$outcomes,
#'                             predictions = PredModelData$predictions)
#' plot(fit)
dca_predictive_model_pair <- function(outcomes,
                                      predictions1,
                                      predictions2,
                                      keep_draws = TRUE,
                                      keep_fit = FALSE,
                                      prior_p = c(1, 1),
                                      prior_se = c(1, 1),
                                      prior_sp = c(1, 1),
                                      summary_probs = c(0.025, 0.975),
                                      refresh = 0, ...) {
  
  df <- get_thr_data(outcomes = outcomes,
                     predictions = predictions)
  f <- function(i, n) lapply(i, function(k) rep(k, n))
  prior_se <- f(prior_se, nrow(df))
  prior_sp <- f(prior_sp, nrow(df))
  
  fit <- .dca_stan(
    N = df$N, d = df$d, tp = df$tp, tn = df$tn, B = nrow(df),
    thresholds = df$thresholds, model_type = "mo",
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    refresh = refresh,
    ...
  )
  
  # TODO: process fit summary for each model
  
  fit_summary <- rstan::summary(
    fit, probs = summary_probs
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("par_name") %>%
    tibble::as_tibble() %>%
    dplyr::select(-se_mean) %>%
    dplyr::rename(estimate := mean)
  
  model_parameters <- fit_summary %>%
    dplyr::filter(
      stringr::str_detect(par_name, "p|Sp|Se")
    )
  net_benefit <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "net_benefit")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = df$thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  treat_all <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "treat_all")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = df$thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  
  output_data <- list(
    model_parameters = model_parameters,
    net_benefit = net_benefit, treat_all = treat_all,
    N = df$N, d = df$d, tp = df$tp, tn = df$tn,
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    thresholds = df$thresholds
  )
  
  if (isTRUE(keep_fit)) {
    output_data[['fit']] <- fit
  }
  
  if(isTRUE(keep_draws)) {
    output_data[['draws']] <- rstan::extract(fit)
  }
  
  .output <- structure(output_data, class = "PredModelDCA")
  return(.output)
}

#' @title Print DiagTestDCA
#'
#' @param obj DiagTestDCA object
#' @importFrom stringr str_glue
#' @export
print.DiagTestDCA <- function(obj, ...) {
  g <- function(obj, .param) {
    info <- obj$model_parameters %>%
      dplyr::filter(par_name == .param)
    paste0(
      .param, ": ",
      round(info$estimate*100, 1), "% [",
      round(info$`2.5%`*100, 1), "% \u2014 ",
      round(info$`97.5%`*100, 1), "%]"
    )
  }
  p <- g(obj = obj, .param = "p")
  se <- g(obj = obj, .param = "Se")
  sp <- g(obj = obj, .param = "Sp")
  .data <- paste0(
    "N=", obj$N, "\tD=", obj$d,
    "\tTP=", obj$tp, "\tTN=", obj$tn
  )
  cat(
    paste0(
      c(
        "DiagTestDCA\n",
        p, se, sp,
        paste0("Number of thresholds: ", length(obj$thresholds)),
        "\nRaw data:", .data
      ),
      collapse = "\n"
    )
  )
}

#' @title Print PredModelDCA
#'
#' @param obj PredModelDCA object
#' @export
print.PredModelDCA <- function(obj, ...) {
  g <- function(obj, .param) {
    info <- obj$model_parameters %>%
      dplyr::filter(par_name == .param)
    paste0(
      .param, ": ",
      round(info$estimate*100, 1), "% [",
      round(info$`2.5%`*100, 1), "% \u2014 ",
      round(info$`97.5%`*100, 1), "%]"
    )
  }
  p <- g(obj = obj, .param = "p")
  .data <- paste0(
    "N=", unique(obj$N), "\tD=", unique(obj$d)
  )
  cat(
    paste0(
      c(
        "PredModelDCA\n",
        p,
        paste0("Number of thresholds: ", length(obj$thresholds)),
        "\nRaw data:", .data
      ),
      collapse = "\n"
    )
  )
}

#' @title Plot DiagTestDCA
#'
#' @param obj DiagTestDCA object
#' @param data_only If TRUE, returns data for ggplot objects (default is FALSE).
#' @param .color Test curve color. Default is "#4DAF4AFF" (green).
#' @importFrom magrittr %>%
#' @export
#' @examples
#' fit <- dca_binary_test(N = 500, d = 83, tp = 77, tn = 378)
#' plot(fit)
plot.DiagTestDCA <- function(obj, data_only = FALSE, .color = "#4DAF4AFF", ...) {
  .colors <- c("Test" = .color,
               "Treat all" = "black", "Treat none" = "gray60")
  
  .p <- obj$net_benefit %>%
    ggplot2::ggplot(ggplot2::aes(thr)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`,
                                      fill = "Test"),
                         alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = estimate, color = "Test")) +
    ggplot2::geom_hline(
      ggplot2::aes(color = "Treat none", yintercept = 0),
      linetype = 'longdash', lwd = 0.8
    ) +
    ggplot2::geom_line(
      data = obj$treat_all,
      ggplot2::aes(thr, estimate, color = "Treat all")
    ) +
    ggplot2::geom_ribbon(
      data = obj$treat_all,
      ggplot2::aes(thr, ymax=`97.5%`, ymin=`2.5%`,
                   fill = "Treat all"),
      alpha = 0.3
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(ylim = c(-0.02, NA)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::scale_color_manual(
      values = .colors
    ) +
    ggplot2::scale_fill_manual(
      values = .colors
    ) +
    ggplot2::labs(x = "Threshold", y = "Net Benefit",
                  color = NULL) +
    ggplot2::guides(fill = 'none')
  
  if (isTRUE(data_only)) return(.p$data)
  
  return(.p)
}

#' @title Plot PredModelDCA
#'
#' @param obj PredModelDCA object
#' @param data_only If TRUE, returns data for ggplot objects (default is FALSE).
#' @param .color Model curve color. Default is "#E41A1CFF" (red).
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(PredModelData)
#' head(PredModelData)
#' fit <- dca_predictive_model(outcomes = PredModelData$outcomes,
#'                             predictions = PredModelData$predictions)
#' plot(fit)
plot.PredModelDCA <- function(obj, data_only = FALSE, .color = "#E41A1CFF", ...) {
  
  .colors <- c("Model" = .color,
               "Treat all" = "black", "Treat none" = "gray60")
  .p <- obj$net_benefit %>%
    ggplot2::ggplot(ggplot2::aes(thr)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`,
                                      fill = "Model"),
                         alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = estimate, color = "Model")) +
    ggplot2::geom_hline(
      ggplot2::aes(color = "Treat none", yintercept = 0),
      linetype = 'longdash', lwd = 0.8
    ) +
    ggplot2::geom_line(
      data = obj$treat_all,
      ggplot2::aes(thr, estimate, color = "Treat all")
    ) +
    ggplot2::geom_ribbon(
      data = obj$treat_all,
      ggplot2::aes(thr, ymax=`97.5%`, ymin=`2.5%`,
                   fill = "Treat all"),
      alpha = 0.3
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(ylim = c(-0.02, NA)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::scale_color_manual(
      values = .colors
    ) +
    ggplot2::scale_fill_manual(
      values = .colors
    ) +
    ggplot2::labs(x = "Threshold", y = "Net Benefit",
                  color = NULL) +
    ggplot2::guides(fill = 'none')
  
  
  if (isTRUE(data_only)) return(.p$data)
  
  return(.p)
}

#' @title Plot DCA list
#'
#' @param ... Fits to plot
#' @param data_only If TRUE, returns data for ggplot objects (default is FALSE).
#' @param .colors Character vector with colors for each fit object.
#' @param .names Optional character vector with names for objects (mostly for internal use).
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(PredModelData)
#' head(PredModelData)
#' # binary tests
#' fit1 <- dca_binary_test(N = 1000, d = 120, tp = 80, tn = 800)
#' fit2 <- dca_binary_test(N = 1000, d = 120, tp = 108, tn = 616)
#' # predictive model
#' data(PredModelData)
#' fit3 <- dca_predictive_model(outcomes = PredModelData$outcomes,
#'                              predictions = PredModelData$predictions)
#' # plot decision curves
#' plot_dca_list("test A" = fit1, "test B" = fit2, "model" = fit3)
#'
#' # names without quotes
#' plot_dca_list(test = fit1, another_test = fit2, model = fit3)
#'
#' # no names (use var names)
#' plot_dca_list(fit1, fit2, fit3)
plot_dca_list <- function(..., data_only = FALSE, .colors = NULL, .names = NULL) {
  fit_list <- list(...)
  if (!is.null(.names)) {
    names(fit_list) <- .names
  } else {
    if (is.null(names(fit_list))) {
      .dots <- match.call(expand.dots = FALSE)$...
      .dots <- .dots[sapply(.dots, is.name)]
      names(fit_list) <- sapply(.dots, deparse)
    }
  }
  
  plot_data <- purrr::map(fit_list, ~ .x$net_benefit) %>%
    dplyr::bind_rows(.id = "fit_name")
  
  ref_fit <- names(fit_list)[1]
  treat_all_data <- fit_list[[ref_fit]]$treat_all %>%
    dplyr::mutate(fit_name = "treat all")
  
  if (is.null(.colors)) {
    n_fits <- length(fit_list)
    .colors <- grDevices:::colorRampPalette(
      paletteer::paletteer_d("RColorBrewer::Dark2", 6)
    )(n_fits)
  }
  names(.colors) <- names(fit_list)
  .colors <- c(.colors, "treat all" = "black", "treat none" = "gray60")
  
  .p <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x = thr)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`,
                                      fill = fit_name),
                         alpha = 0.3) +
    ggplot2::geom_ribbon(data = treat_all_data,
                         ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`),
                         alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = estimate, color = fit_name)) +
    ggplot2::geom_line(data = treat_all_data,
                       ggplot2::aes(y = estimate, color = fit_name)) +
    ggplot2::geom_hline(
      ggplot2::aes(color = "treat none", yintercept = 0),
      linetype = 'longdash', lwd = 0.8
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(ylim = c(-0.02, NA)) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::scale_color_manual(
      values = .colors,
      breaks = names(.colors)
    ) +
    ggplot2::scale_fill_manual(
      values = .colors,
      breaks = names(.colors)
    ) +
    ggplot2::labs(x = "Threshold", y = "Net Benefit",
                  color = NULL)
  
  if (isTRUE(data_only)) return(.p$data)
  
  return(.p)
}

#' @title Get Threshold Performance Data
#'
#' @param outcomes Integer vector (0 or 1) with binary outcomes.
#' @param predictions Numeric vector with predicted probabilities.
#' @importFrom magrittr %>%
get_thr_data <- function(outcomes,
                         predictions,
                         thresholds = seq(0, 0.5, 0.02)) {
  
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

#' @title Plot Delta Net Benefit
#' @importFrom magrittr %>%
#' @import patchwork
example_delta_nb <- function() {
  
  as_list = FALSE
  
  d1 <- simulate_diagnostic_test_data(B = 1, N = 5000,
                                      true_p = 0.1,
                                      true_se = 0.95,
                                      true_sp=0.7)
  d2 <- simulate_diagnostic_test_data(B = 1, N = 5000,
                                      true_p = 0.1,
                                      true_se = 0.9,
                                      true_sp=0.8)
  fit1 <- dca_binary_test(N = d1$N,
                          d = d1$d,
                          tp = d1$tp,
                          tn = d1$tn)
  fit2 <- dca_binary_test(N = d2$N,
                          d = d2$d,
                          tp = d2$tp,
                          tn = d2$tn)
  
  
  true_nb1 <- sapply(fit1$thresholds,function(i) {
    .1*.95 - (1-.1)*(1-.7)*(i/(1-i))
  })
  true_nb2 <- sapply(fit2$thresholds,function(i) {
    .1*.9 - (1-.1)*(1-.8)*(i/(1-i))
  })
  true_delta <- purrr::map2_dbl(true_nb2, true_nb1, ~{
    .x-.y
  })
  
  d <- fit2$draws$net_benefit - fit1$draws$net_benefit
  q <- matrixStats::colQuantiles(d, probs = c(.025, .975))
  dp <- tibble::tibble(
    thr = fit1$thresholds,
    estimate = matrixStats::colMedians(d),
    `2.5%` = q[,'2.5%'],
    `97.5%` = q[,'97.5%']
  ) %>%
    ggplot2::ggplot(ggplot2::aes(thr, y=estimate,ymin=`2.5%`, ymax=`97.5%`)) +
    ggplot2::geom_ribbon(alpha=.3) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = true_delta), col='red'
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 'longdash',
      color = 'gray30', lwd = 0.8
    ) +
    ggplot2::labs(x = "Threshold", y = "NB2 - NB1")
  
  prob_better <- tibble::tibble(
    thr = fit1$thresholds,
    estimate = colMeans(d > 0)
  ) %>%
    ggplot2::ggplot(ggplot2::aes(thr, y=estimate)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(1),
      breaks = scales::pretty_breaks(10)
    ) +
    ggplot2::labs(x = "Threshold", y = "Pr(NB2 > NB1)")
  
  dca <- bayesDCA::plot_dca_list(test1 = fit1, test2 = fit2)
  
  if (isTRUE(as_list)) {
    return(list(dca = dca, delta = dp))
  }
  
  dca | dp
}


#' @title Compare DCA
#' @param ... Pass two bayesDCA fit objects with names as they should appear in plots'
#' legends (e.g. "Test one" = fit1, test2 = fit2).
#' @param as_plot_list If TRUE, returns list of ggplot2 plots (default is FALSE).
#' @param .names Optional character vector with names for objects (mostly for internal use).
#' @return A figure of multiple plots
#' ([`patchwork`](https://patchwork.data-imaginist.com/) object).
#' @examples
#'
#'# Test 1: more specificity
#'data1 <- simulate_diagnostic_test_data(B = 1, N = 5000,
#'                                       true_p = 0.1, true_se = 0.8, true_sp = 0.8)
#'head(data1)
#'
#'# Test 2: more sensitivity
#'data2 <- simulate_diagnostic_test_data(B = 1, N = 5000,
#'                                       true_p = 0.1, true_se = 0.95, true_sp = 0.7)
#'
#'fit1 <- dca_binary_test(N = data1$N, d = data1$d, tp = data1$tp, tn = data1$tn)
#'fit2 <- dca_binary_test(N = data2$N, d = data2$d, tp = data2$tp, tn = data2$tn)
#'compare_dca("More Spec" = fit1, "More Sens" = fit2)
#'
#' @importFrom magrittr %>%
#' @import patchwork
#' @export

compare_dca <- function(..., as_plot_list = FALSE, .names = NULL) {
  dots <- list(...)
  if (!is.null(.names)) {
    names(dots) <- .names
  } else {
    if (is.null(names(dots))) {
      .dots <- match.call(expand.dots = FALSE)$...
      .dots <- .dots[sapply(.dots, is.name)]
      names(dots) <- sapply(.dots, deparse)
    }
  }
  
  dp <- plot_delta_nb(..., .names = names(dots))
  prob_better <- plot_prob_better(..., .names = names(dots))
  dca <- plot_dca_list(..., .names = names(dots))
  
  if (isTRUE(as_plot_list)) {
    return(list(dca = dca, delta = dp, prob_better = prob_better))
  }
  
  dca / (dp | prob_better)
}


#' @title Plot delta net benefit
#' @param ... Pass two bayesDCA fit objects with names as they should appear in plots'
#' legends (e.g. "Test one" = fit1, test2 = fit2).
#' @param data_only If TRUE, returns data for ggplot objects (default is FALSE).
#' @param .delta Optional pre-computed delta net benefit (mostly for internal use).
#' @param .names Optional character vector with names for objects (mostly for internal use).
#' @return A ggplot object.

plot_delta_nb <- function(..., data_only = FALSE, .delta = NULL, .names = NULL) {
  dots <- list(...)
  n_dots <- length(dots)
  if (!is.null(.names)) {
    names(dots) <- .names
  } else {
    if (is.null(names(dots))) {
      .dots <- match.call(expand.dots = FALSE)$...
      .dots <- .dots[sapply(.dots, is.name)]
      names(dots) <- sapply(.dots, deparse)
    }
  }
  
  if (is.null(.delta)) {
    if (n_dots == 1) {
      .delta <- dots[[1]]$draws$delta
    } else if (n_dots == 2) {
      .delta <- dots[[1]]$draws$net_benefit - dots[[2]]$draws$net_benefit
    } else {
      stop("Cannot plot more than two objects.")
    }
    
  }
  
  q <- matrixStats::colQuantiles(.delta, probs = c(.025, .5, .975))
  
  .subtitle <- ifelse(
    n_dots == 1,
    paste0(names(dots)[1], ' \u2212 Treat all/none'),
    paste0(names(dots)[1], ' \u2212 ', names(dots)[2])
  )
  
  .plot <- tibble::tibble(
    thr = dots[[1]]$thresholds,
    estimate = q[,'50%'],
    `2.5%` = q[,'2.5%'],
    `97.5%` = q[,'97.5%']
  ) %>%
    ggplot2::ggplot(ggplot2::aes(thr, y=estimate, ymin=`2.5%`, ymax=`97.5%`)) +
    ggplot2::geom_ribbon(alpha=.3) +
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
      x = "Threshold",
      y = '\u0394 NB',
      subtitle = .subtitle
    )
  
  if (isTRUE(data_only)) return(.plot$data)
  
  return(.plot)
}

#' @title Plot probability of better net benefit between two models or tests
#' @param ... Pass two bayesDCA fit objects with names as they should appear in plots'
#' legends (e.g. "Test one" = fit1, test2 = fit2).
#' @param data_only If TRUE, returns data for ggplot objects (default is FALSE).
#' @param .delta Optional pre-computed delta net benefit (mostly for internal use).
#' @param .names Optional character vector with names for objects (mostly for internal use).
#' @return A ggplot object.
plot_prob_better <- function(..., data_only = FALSE, .delta = NULL, .names = NULL) {
  dots <- list(...)
  n_dots <- length(dots)
  if (!is.null(.names)) {
    names(dots) <- .names
  } else {
    if (is.null(names(dots))) {
      .dots <- match.call(expand.dots = FALSE)$...
      .dots <- .dots[sapply(.dots, is.name)]
      names(dots) <- sapply(.dots, deparse)
    }
  }
  
  if (is.null(.delta)) {
    if (n_dots == 1) {
      .delta <- dots[[1]]$draws$delta
    } else if (n_dots == 2) {
      .delta <- dots[[1]]$draws$net_benefit - dots[[2]]$draws$net_benefit
    } else {
      stop("Cannot plot more than two objects.")
    }
  }
  
  .subtitle <- ifelse(
    n_dots == 1,
    paste0("Pr( ", names(dots)[1], " > Treat all/none", " )"),
    paste0("Pr( ", names(dots)[1], " > ", names(dots)[2], " )")
  )
  
  .plot <- tibble::tibble(
    thr = dots[[1]]$thresholds,
    estimate = colMeans(.delta > 0)
  ) %>%
    ggplot2::ggplot(ggplot2::aes(thr, y = estimate)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(1),
      breaks = scales::pretty_breaks(10)
    ) +
    ggplot2::labs(
      x = "Threshold", y = NULL,
      subtitle = .subtitle
    )
  
  if (isTRUE(data_only)) return(.plot$data)
  
  return(.plot)
}

#' @title Compare posterior distributions at specific thresholds
#' @param ... Pass up to two bayesDCA fit objects with names as they should appear in plots'
#' legends (e.g. "Test one" = fit1, test2 = fit2). First object appears in y axis.
#' @param thr Threshold to use - should have been computed in the first object.
#' @param .names Optional character vector with names for objects (mostly for internal use).
#' @param data_only If TRUE, returns data for ggplot objects (default is FALSE).
#' @return A ggplot object.
compare_threshold_posterior <- function(..., thr, .names = NULL, data_only = FALSE) {
  dots <- list(...)
  n_dots <- length(dots)
  if (!is.null(.names)) {
    names(dots) <- .names
  } else {
    if (is.null(names(dots))) {
      .dots <- match.call(expand.dots = FALSE)$...
      .dots <- .dots[sapply(.dots, is.name)]
      names(dots) <- sapply(.dots, deparse)
    }
  }
  thr_ix <- which(dplyr::near(dots[[1]]$thresholds, thr))[1]
  if (is.na(thr_ix)) {
    stop("There seems to have no such threshold available.")
  }
  
  yvalues <- dots[[1]]$draws$net_benefit[, thr_ix]
  if (n_dots == 1) {
    xvalues <- dots[[1]]$draws$treat_all[, thr_ix]
    xlab <- "Treat all"
  } else if (n_dots == 2) {
    xvalues <- dots[[2]]$draws$net_benefit[, thr_ix]
    xlab <- names(dots)[2]
  } else {
    stop("Cannot plot more than two objects.")
  }
  
  diff_prob <- round(mean(yvalues > xvalues)*100, 2)
  .subtitle <- paste0("Pr( ", names(dots)[1], " > ", xlab, " ) = ", diff_prob, "%")
  .plot <- data.frame(xvalues, yvalues) %>%
    ggplot2::ggplot(ggplot2::aes(xvalues, yvalues)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         color = "red", linetype = "longdash") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = xlab,
      y = names(dots)[1],
      title = paste0("Threshold: ", thr),
      subtitle = .subtitle
    )
  
  if (isTRUE(data_only)) return(.plot$data)
  
  return(.plot)
}
