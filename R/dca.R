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
                      model_type = "mo",
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
                                model_type = "correlated",
                                prior_p = c(1, 1),
                                prior_se = c(1, 1),
                                prior_sp = c(1, 1),
                                refresh = 0, ...) {

  fit <- .dca_stan(
    N = N, d = d, tp = tp, tn = tn, model_type = model_type,
    prior_p = prior_p, prior_se = prior_se, prior_sp = prior_sp,
    thresholds = thresholds, refresh = refresh,
    ...
  )

  fit_summary <- rstan::summary(
    fit, probs = c(0.025, 0.1, .25, .75, 0.9, 0.975)
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("par_name") %>%
    tibble::as_tibble() %>%
    dplyr::select(-se_mean)

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
    fit, probs = c(0.025, 0.1, .25, .75, 0.9, 0.975)
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("par_name") %>%
    tibble::as_tibble() %>%
    dplyr::select(-se_mean)

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
      round(info$mean*100, 1), "% [",
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
