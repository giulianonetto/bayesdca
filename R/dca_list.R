#' Bayesian Decision Curve Analysis for a list of Predictive Models
#'
#' @export
#' @description Estimate decision curves for a list of predictive models
#' all at once. This is necessary to make comparative inferences across models
#' using their corresponding posterior draws.
#' @param .data A data.frame with an `outcomes` column (0 or 1 for each individual)
#' and one or more columns with predicted probabilities from each of desired list
#' of predictive models.
#' @param thresholds Numeric vector with probability thresholds with which
#' the net benefit should be computed (default is `seq(0, 0.5, 0.01)`).
#' @param keep_fit Logical indicating whether to keep `stanfit` in
#' the output (default is FALSE).
#' @param keep_draws Logical indicating whether to keep posterior
#' draws from `stanfit` object (default is TRUE).
#' @param prior_p,prior_se,prior_sp Non-negative shape values for
#' Beta(alpha, beta) priors used for p, Se, and Sp, respectively. 
#' Default is uniform prior for all parameters - Beta(1, 1).
#' For `prior_p`, a single vector of the form `c(a, b)` can be provided.
#' For both `prior_se` and `prior_sp`, the prior can be specified for each
#' model and is internally propagated through all thresholds. Input type
#' should be a 2-column matrix with one row per input model, 
#' in the same order as the prediction columns of `.data`. Each row represents `c(a, b)`
#' for each corresponding model - first row for first model and so on.
#' @param refresh Control verbosity of `rstan::sampling` (check its help
#' page for details).
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `PredModelListDCA`
#' @importFrom magrittr %>%
#' @examples
#' data(PredModelData)
#' head(PredModelData)
#' fit <- dca_predictive_model(.data = PredModelData)
#' plot(fit)
dca_predictive_model <- function(.data,
                                 thresholds = seq(0, 0.5, 0.01),
                                 keep_draws = TRUE,
                                 keep_fit = FALSE,
                                 prior_p = NULL,
                                 prior_se = NULL,
                                 prior_sp = NULL,
                                 summary_probs = c(0.025, 0.975),
                                 refresh = 0, 
                                 ...) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }
  model_names <- colnames(.data)[-1]
  threshold_data <- get_thr_data_list(.data = .data,
                                      thresholds = thresholds)
  priors <- get_prior_parameters(prior_p = prior_p,
                                 prior_se = prior_se,
                                 prior_sp = prior_sp)
  
  fit <- .dca_stan_list(
    n_thr = length(thresholds),
    n_models = ncol(.data) - 1,
    N = threshold_data$N,
    d = threshold_data$d,
    tp = threshold_data$tp,
    tn = threshold_data$tn,
    thresholds = threshold_data$thresholds,
    prior_p1 = priors[['p1']],
    prior_p2 = priors[['p2']],
    prior_Se1 = priors[['Se1']],
    prior_Se2 = priors[['Se2']],
    prior_Sp1 = priors[['Sp1']],
    prior_Sp2 = priors[['Sp2']],
    refresh = refresh,
    ...
  )
  
  dca_list <- extract_dca_list(fit = fit,
                               summary_probs = summary_probs,
                               thresholds = thresholds,
                               model_names = model_names,
                               keep_draws = keep_draws)
  
  output_data <- list(
    dca_list = dca_list,
    thresholds = thresholds,
    delta = delta,
    .data = .data,
    threshold_data = threshold_data,
    priors = priors,
    model_names = model_names,
  )
  
  if (isTRUE(keep_fit)) {
    output_data[['fit']] <- fit
  }
  
  .output <- structure(output_data, class = "PredModelListDCA")
  return(.output)
}

get_prior_parameters <- function() {
  f <- function(i, n) lapply(i, function(k) rep(k, n))
  prior_se <- f(prior_se, nrow(df))
  prior_sp <- f(prior_sp, nrow(df))
}

extract_fit_summary <- function(fit,
                                summary_probs) {
  fit_summary <- rstan::summary(
    fit, probs = summary_probs
  )$summary %>%
    data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("par_name") %>%
    tibble::as_tibble() %>%
    dplyr::select(-se_mean) %>%
    dplyr::rename(estimate := mean)
  
  return(fit_summary)
}

extract_model_parameters <- function(fit_summary) {
  
  model_parameters <- fit_summary %>%
    dplyr::filter(
      stringr::str_detect(par_name, "p|Sp|Se")
    )
  return(model_parameters)
  
}

extract_net_benefit <- function(fit_summary,
                                thresholds) {
  net_benefit <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "net_benefit")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
  return(net_benefit)
  
}



get_fit_summaries <- function(fit, summary_probs, thresholds) {
  
  
  
  treat_all <- fit_summary %>%
    dplyr::filter(stringr::str_detect(par_name, "treat_all")) %>%
    dplyr::mutate(
      i = as.numeric(stringr::str_extract(par_name, "\\d+")),
      thr = thresholds[i]
    ) %>%
    dplyr::select(par_name, thr, dplyr::everything(), -i)
}