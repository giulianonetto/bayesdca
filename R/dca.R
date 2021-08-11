#' Fit Bayesian Decision Curve Analysis using Stan
#'
#' @export
#' @param N Integer specifying total number of samples (i.e., participants).
#' @param d Diseased: integer specifying total number of diseased persons.
#' @param tp True Positives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model.
#' @param tn True Negatives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param model_type Character indicating model type. Either "correlated"
#' (default) or "independent".
#' @param prior_p,prior_se,prior_sp Numeric vector with shapes for
#' Beta(alpha, beta) prior used for p, Se, and Sp. Default is c(1, 1) for all.
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
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

#' Bayesian Decision Curve Analysis for Diagnostic Test
#'
#' @export
#' @param N Integer specifying total number of samples (i.e., participants).
#' @param d Diseased: integer specifying total number of diseased persons.
#' @param tp True Positives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model.
#' @param tn True Negatives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default FALSE).
#' @param model_type Character indicating model type. Either "correlated"
#' (default) or "independent".
#' @param prior_p,prior_se,prior_sp Numeric vector with shapes for
#' Beta(alpha, beta) prior used for p, Se, and Sp. Default is c(1, 1) for all.
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `DiagTestDCA`
#' @importFrom magrittr %>%
#'
dca_diagnostic_test <- function(N, d, tp, tn,
                                thresholds = seq(0, 0.5, 0.01),
                                keep_fit = FALSE,
                                keep_draws = TRUE,
                                model_type = "ic",
                                prior_p = c(1, 1),
                                prior_se = c(1, 1),
                                prior_sp = c(1, 1),
                                summary_probs = c(0.025, 0.975),
                                refresh = 0, ...) {

  thresholds <- pmin(thresholds, 0.99)

  fit <- .dca_stan(
    N = N, d = d, tp = tp, tn = tn, model_type = model_type,
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
    thresholds = thresholds,
    model_type = model_type
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
#' @param N Integer specifying total number of samples (i.e., participants).
#' @param d Diseased: integer specifying total number of diseased persons.
#' @param tp True Positives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model.
#' @param tn True Negatives: integer specifying total number of diseased persons correctly
#' identified as such by the diagnostic test of prediction model.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default FALSE).
#' @param model_type Character indicating model type. Either "correlated"
#' (default) or "independent".
#' @param prior_p,prior_se,prior_sp Numeric vector with shapes for
#' Beta(alpha, beta) prior used for p, Se, and Sp. Default is c(1, 1) for all.
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `DiagTestDCA`
#' @importFrom magrittr %>%
#'
dca_predictive_model <- function(df,
                                 keep_fit = FALSE,
                                 model_type = "mo",
                                 prior_p = c(1, 1),
                                 prior_se = c(1, 1),
                                 prior_sp = c(1, 1),
                                 summary_probs = c(0.025, 0.975),
                                 refresh = 0, ...) {

  f <- function(i, n) lapply(i, function(k) rep(k, n))
  prior_se <- f(prior_se, nrow(df))
  prior_sp <- f(prior_sp, nrow(df))

  fit <- .dca_stan(
    N = df$N, d = df$d, tp = df$tp, tn = df$tn, B = nrow(df),
    thresholds = df$thresholds, model_type = model_type,
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
    thresholds = df$thresholds,
    model_type = model_type
  )

  if (isTRUE(keep_fit)) {
    output_data[['fit']] <- fit
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
        paste0("DiagTestDCA", " ~ ", obj$model_type, " parameters model\n"),
        p, se, sp,
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
#' @param type Either "net benefit" (alias "nb" or "dca") or "parameters"
#' (alias "p"). If net benefit, plots estimated net benefit against
#' thresholds (DCA); if parameters, plots estimated prevalence,
#' sensitivity, and specificity.
#' @importFrom magrittr %>%
#' @export
plot.DiagTestDCA <- function(obj, type = "nb", just_df=FALSE, ...) {
  if (type %in% c("nb", "net benefit")) {
    .p <- obj$net_benefit %>%
      ggplot2::ggplot(ggplot2::aes(thr)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`),
                           alpha = .5) +
      ggplot2::geom_line(ggplot2::aes(y = estimate)) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = 'longdash',
        color = 'gray30', lwd = 0.8
      ) +
      ggplot2::geom_line(
        data = obj$treat_all,
        ggplot2::aes(thr, estimate, color = "Treat all")
      ) +
      ggplot2::geom_ribbon(
        data = obj$treat_all, alpha = .3, fill = "red",
        ggplot2::aes(thr, ymax=`97.5%`, ymin=`2.5%`)
      ) +
      ggplot2::theme_bw() +
      ggplot2::coord_cartesian(ylim = c(-0.02, NA)) +
      ggplot2::scale_x_continuous(
        labels = scales::percent_format(1)
      ) +
      ggplot2::labs(x = "Threshold", y = "Net Benefit",
                    color = NULL)
  } else if (type %in% c("p", "parameters")) {
    .p <- obj$model_parameters %>%
      ggplot2::ggplot(ggplot2::aes(x = par_name)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(
          y = estimate,
          ymin = `2.5%`, ymax = `97.5%`
        )
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::scale_y_continuous(labels = scales::percent,
                                  breaks = scales::pretty_breaks()) +
      ggplot2::labs(x = NULL, y = NULL)
  } else {
    stop("Invalid 'type' argument.")
  }

  if (isTRUE(just_df)) return(.p$data)

  return(.p)
}

#' @title Plot PredModelDCA
#'
#' @param obj PredModelDCA object
#' @param type Either "net benefit" (alias "nb" or "dca") or "parameters"
#' (alias "p"). If net benefit, plots estimated net benefit against
#' thresholds (DCA); if parameters, plots estimated prevalence,
#' sensitivity, and specificity.
#' @importFrom magrittr %>%
#' @export
plot.PredModelDCA <- function(obj, type = "nb", just_df=FALSE, ...) {
  if (type %in% c("nb", "net benefit")) {
    .p <- obj$net_benefit %>%
      ggplot2::ggplot(ggplot2::aes(thr)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`),
                           alpha = .5) +
      ggplot2::geom_line(ggplot2::aes(y = mean)) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = 'longdash',
        color = 'gray30', lwd = 0.8
      ) +
      ggplot2::geom_line(
        data = obj$treat_all,
        ggplot2::aes(thr, mean, color = "Treat all")
      ) +
      ggplot2::geom_ribbon(
        data = obj$treat_all, alpha = .3, fill = "red",
        ggplot2::aes(thr, ymax=`97.5%`, ymin=`2.5%`)
      ) +
      ggplot2::theme_bw() +
      ggplot2::coord_cartesian(ylim = c(-0.02, NA)) +
      ggplot2::scale_x_continuous(
        labels = scales::percent_format(1)
      ) +
      ggplot2::scale_y_continuous(
        breaks = scales::pretty_breaks()
      ) +
      ggplot2::labs(x = "Threshold", y = "Net Benefit",
                    color = NULL)
  } else if (type %in% c("p", "parameters")) {
    .p <- obj$model_parameters %>%
      ggplot2::ggplot(ggplot2::aes(x = par_name)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(
          y = mean,
          ymin = `2.5%`, ymax = `97.5%`
        )
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::scale_y_continuous(labels = scales::percent,
                                  breaks = scales::pretty_breaks()) +
      ggplot2::labs(x = NULL, y = NULL)
  } else {
    stop("Invalid 'type' argument.")
  }

  if (isTRUE(just_df)) return(.p$data)

  return(.p)
}

#' @title Plot DCA list
#'
#' @param ... Fits to plot
#' @param type Either "net benefit" (alias "nb" or "dca") or "parameters"
#' (alias "p"). If net benefit, plots estimated net benefit against
#' thresholds (DCA); if parameters, plots estimated prevalence,
#' sensitivity, and specificity.
#' @importFrom magrittr %>%
#' @export
plot_dca_list <- function(..., type = "nb", data_only=FALSE) {
  if (type %in% c("nb", "net benefit")) {
    fit_list <- list(...)

    if (is.null(names(fit_list))) {
      stop("Fits must have names.")
    }

    plot_data <- purrr::map(fit_list, ~ .x$net_benefit) %>%
      dplyr::bind_rows(.id = "fit_name")

    ref_fit <- names(fit_list)[1]
    treat_all_data <- fit_list[[ref_fit]]$treat_all %>%
      dplyr::mutate(fit_name = "treat all")

    .colors <- paletteer::paletteer_d("RColorBrewer::Dark2", n = length(fit_list))
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
  }  else {
    stop("Invalid 'type' argument.")
  }

  if (isTRUE(data_only)) return(.p$data)

  return(.p)
}

#' @title Get Threshold Performance Data
#'
#' @param df Input data.frame
#' @param outcome Character indicating the column with binary outcome (as 0/1s).
#' @param prediction Character indicating the column with predicted probabilities.
#' @importFrom magrittr %>%
#' @export
get_thr_data <- function(df,
                         outcome = "outcomes",
                         prediction = "predictions",
                         thresholds = seq(0, 0.5, 0.02)) {

  thr_data <- tibble::tibble(
    N = nrow(df),
    d = sum(df[[outcome]]),
    thresholds = thresholds
  ) %>%
    dplyr:::mutate(
      thr_perf = purrr::map(thresholds, function(.thr) {
        tp <- sum(df[df[[outcome]]==1,][[prediction]] >= .thr)
        tn <- sum(df[df[[outcome]]==0,][[prediction]] < .thr)
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
#' @export
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
  fit1 <- dca_diagnostic_test(N = d1$N,
                              d = d1$d,
                              tp = d1$tp,
                              tn = d1$tn)
  fit2 <- dca_diagnostic_test(N = d2$N,
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
    estimate = colMeans(d),
    `2.5%` = q[,'2.5%'],
    `97.5%` = q[,'97.5%']
  ) %>%
    ggplot(aes(thr, y=estimate,ymin=`2.5%`, ymax=`97.5%`)) +
    geom_ribbon(alpha=.3) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    geom_line(
      aes(y = true_delta), col='red'
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
    ggplot(aes(thr, y=estimate)) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    scale_y_continuous(
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


#' @title Compare Net Benefit
#' @param ... Pass two bayesDCA fit objects with names as they should appear in plots'
#' legends (e.g. "Test one" = fit1, test2 = fit2).
#' @param as_plot_list If TRUE, returns list of ggplot2 plots (default is FALSE).
#' @return A figure of multiple plots
#' ([`patchwork`](https://patchwork.data-imaginist.com/) object).
#' @examples
#'
#'# Test 1: more specificity
#'data1 <- simulate_diagnostic_test_data(B = 1, N = 5000,
#'                                       true_p = 0.1, true_se = 0.9, true_sp = 0.8)
#'# Test 2: more sensitivity
#'data2 <- simulate_diagnostic_test_data(B = 1, N = 5000,
#'                                       true_p = 0.1, true_se = 0.95, true_sp = 0.7)
#'
#'fit1 <- dca_diagnostic_test(N = data1$N, d = data1$d, tp = data1$tp, tn = data1$tn)
#'fit2 <- dca_diagnostic_test(N = data2$N, d = data2$d, tp = data2$tp, tn = data2$tn)
#'compare_net_benefit("More Spec" = fit1, "More Sens" = fit2)
#'
#' @importFrom magrittr %>%
#' @import patchwork
#' @export

compare_net_benefit <- function(..., as_plot_list = FALSE) {

  dots <- list(...)
  d <- dots[[1]]$draws$net_benefit - dots[[2]]$draws$net_benefit
  q <- matrixStats::colQuantiles(d, probs = c(.025, .975))
  dp <- tibble::tibble(
    thr = dots[[1]]$thresholds,
    estimate = colMeans(d),
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
    ggplot2::geom_hline(
      yintercept = 0, linetype = 'longdash',
      color = 'gray30', lwd = 0.8
    ) +
    ggplot2::labs(
      x = "Threshold",
      y = '\u0394 NB',
      subtitle = paste0(names(dots)[1], ' \u2212 ', names(dots)[2])
    )

  prob_better <- tibble::tibble(
    thr = dots[[1]]$thresholds,
    estimate = colMeans(d > 0)
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
      subtitle = paste0("Pr( ", names(dots)[1], " > ", names(dots)[2], " )")
    )

  dca <- bayesDCA::plot_dca_list(...)

  if (isTRUE(as_plot_list)) {
    return(list(dca = dca, delta = dp, prob_better = prob_better))
  }

  dca / (dp | prob_better)
}
