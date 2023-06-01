#' @title Fit Bayesian Decision Curve Analysis
#' using Stan for survival outcomes (Weibull model)
#'
#' @param refresh Control verbosity of
#' [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).  # nolint
#' @param ... Arguments passed to
#' [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).  # nolint
#' @return An object of class
#' [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)  # nolint
#' @keywords internal
.dca_stan_surv_weibull <- function(sample_size, # nolint
                                   n_thr,
                                   n_strategies,
                                   event_times_stacked,
                                   censored_times_stacked,
                                   n_event_times_stacked,
                                   n_censored_times_stacked,
                                   event_times_start_positions,
                                   censored_times_start_positions,
                                   event_times_sizes,
                                   censored_times_sizes,
                                   prediction_time,
                                   thresholds,
                                   pos_post1,
                                   pos_post2,
                                   total_event_times,
                                   total_censored_times,
                                   event_times_marginal,
                                   censored_times_marginal,
                                   other_models_indices,
                                   prior_only,
                                   shape_prior,
                                   scale_prior,
                                   shape_prior_pars,
                                   scale_prior_pars,
                                   iter = 2000,
                                   refresh = 0,
                                   ...) {
  thresholds <- pmin(thresholds, 0.9999) # nolint odds(1) = Inf

  standata <- list(
    N = sample_size,
    n_thr = n_thr,
    n_strategies = n_strategies,
    event_times_stacked = event_times_stacked,
    censored_times_stacked = censored_times_stacked,
    n_event_times_stacked = n_event_times_stacked,
    n_censored_times_stacked = n_censored_times_stacked,
    event_times_start_positions = event_times_start_positions,
    censored_times_start_positions = censored_times_start_positions,
    event_times_sizes = event_times_sizes,
    censored_times_sizes = censored_times_sizes,
    prediction_time = prediction_time,
    thresholds = thresholds,
    pos_post1 = pos_post1,
    pos_post2 = pos_post2,
    total_event_times = total_event_times,
    total_censored_times = total_censored_times,
    event_times_marginal = event_times_marginal,
    censored_times_marginal = censored_times_marginal,
    other_models_indices = other_models_indices,
    prior_only = prior_only,
    shape_prior = shape_prior,
    scale_prior = scale_prior,
    shape_prior_pars = shape_prior_pars,
    scale_prior_pars = scale_prior_pars,
    iter = iter,
    refresh = refresh
  )

  dots <- list(...)
  if ("control" %in% names(dots)) {
    control <- dots[["control"]]
  } else {
    control <- list(adapt_delta = 0.9)
  }

  .model <- stanmodels$dca_survival_weibull
  stanfit <- rstan::sampling(
    .model,
    data = standata,
    control = control,
    iter = iter,
    refresh = refresh, ...
  )
  return(stanfit)
}

#' @title Bayesian Decision Curve Analysis
#' for Survival outcomes
#' @param .data dataframe whose first column named "outcomes" is a `survival::Surv` object
#' and remaining columns are the decision strategies to assess.
#' @param prediction_time Prediction time horizon (e.g., if models predict risk
#' of death at one year and data is in year, `prediction_time` should be `1`.)
#' @param thresholds Decision thresholds -- within interval (0, 1).
#' @param keep_draws If true, posterior draws are kept in the output object.
#' @param keep_fit If true, `stanfit` object is kept in the output object.
#' @param summary_probs Probabilities for posterior credible intervals (defaults to a 95% Cr.I.).
#' @param positivity_prior Shape parameters for prior on positivity probability.
#' @param shape_prior type of prior distribution for shape parameter of the Weibull distribution. Either "student" or "gamma".
#' @param scale_prior type of prior distribution for scale parameter of the Weibull distribution. Either "student" or "gamma".
#' @param shape_prior_pars vector with prior parameters for the prior shape of the Weibull distribution.
#' If `shape_prior="student"`, it should be a vector of length 3 with degrees of freedom, mean, and scale,
#' respectively; if `shape_prior="gamma"`, it should be a vector of length 2 with shape and rate, respectively.
#' @param prior_only If TRUE, samples from the prior only.
#' @param iter Passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html). Number of iterations/draws for Stan.
#' @param refresh Controls verbosity of
#' [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
#' @param ... Arguments passed to
#' [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).  # nolint
#' @export
#' @return An object of class `BayesDCASurv`
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' data(dca_survival_data)
#' fit <- dca_surv(dca_survival_data, prediction_time = 1, iter = 1000, chains = 1)
#' plot(fit)
#' }
dca_surv <- function(.data, # nolint
                     prediction_time,
                     thresholds = seq(0, 0.5, length = 51),
                     keep_draws = TRUE,
                     keep_fit = FALSE,
                     summary_probs = c(0.025, 0.975),
                     positivity_prior = c(1, 1),
                     shape_prior = c("student", "gamma"),
                     scale_prior = c("student", "gamma"),
                     shape_prior_pars = c(10, 0, 1.5),
                     scale_prior_pars = c(30, 0, 100),
                     prior_only = FALSE,
                     iter = 4000,
                     refresh = 0,
                     ...) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }

  stopifnot(
    "'outcomes' column must be a Surv object. " = survival::is.Surv(.data[["outcomes"]]) # nolint
  )
  stopifnot(
    "positivity_prior must be a vector of length 2." = length(positivity_prior) == 2L # nolint
  )

  shape_prior <- ifelse(
    match.arg(shape_prior) == "student", 1, 2
  )
  scale_prior <- ifelse(
    match.arg(scale_prior) == "student", 1, 2
  )
  shape_prior_pars <- c(shape_prior_pars, 0)[1:3] # fix to size 3
  scale_prior_pars <- c(scale_prior_pars, 0)[1:3]

  # avoid thresholds in {0, 1}
  thresholds <- thresholds %>%
    pmin(0.99) %>%
    pmax(1e-9) %>%
    unique()

  strategies <- colnames(.data)[-1]
  prediction_data <- data.frame(.data[, -1])
  colnames(prediction_data) <- strategies
  n_strategies <- length(strategies)
  n_thresholds <- length(thresholds)
  surv_data <- data.frame(
    .time = unname(.data[["outcomes"]][, 1]), # observed time-to-event
    .status = unname(.data[["outcomes"]][, 2]) # 1 if event, 0 if censored
  )

  event_times_marginal <- surv_data$.time[surv_data$.status == 1L] # nolint
  censored_times_marginal <- surv_data$.time[surv_data$.status == 0L] # nolint
  total_event_times <- length(event_times_marginal)
  total_censored_times <- length(censored_times_marginal)

  event_times_stacked <-
    censored_times_stacked <-
    numeric() # hard to guess the size

  # these are used within Stan with the segment() function
  event_times_sizes <-
    censored_times_sizes <- matrix(
      nrow = n_strategies, ncol = n_thresholds
    )

  event_times_start_positions <-
    censored_times_start_positions <- matrix(
      nrow = n_strategies, ncol = n_thresholds
    )

  event_times_position <- censored_times_position <- 1
  for (j in seq_len(n_strategies)) {
    for (m in seq_len(n_thresholds)) {
      .thr <- thresholds[m]
      .model <- strategies[j]
      .preds <- prediction_data[, .model]

      .event_times_positives <- surv_data$.time[
        .preds >= .thr & surv_data$.status == 1L
      ]
      .censored_times_positives <- surv_data$.time[
        .preds >= .thr & surv_data$.status == 0L
      ]
      .event_times_size <- length(.event_times_positives)
      .censored_times_size <- length(.censored_times_positives)
      event_times_stacked <- c(
        event_times_stacked, .event_times_positives
      )
      censored_times_stacked <- c(
        censored_times_stacked, .censored_times_positives
      )

      event_times_start_positions[j, m] <- event_times_position
      censored_times_start_positions[j, m] <- censored_times_position
      event_times_sizes[j, m] <- .event_times_size
      censored_times_sizes[j, m] <- .censored_times_size

      event_times_position <- event_times_position + .event_times_size
      censored_times_position <- censored_times_position + .censored_times_size
    }
  }

  posterior_positivity_pars <- get_positivity_posterior_parameters(
    .prediction_data = prediction_data,
    .thresholds = thresholds,
    .prior_shape1 = positivity_prior[1],
    .prior_shape2 = positivity_prior[2],
    .prior_only = prior_only
  )

  # random data information
  other_models_indices <- lapply(
    1:n_strategies,
    function(i) (1:n_strategies)[-i]
  )

  fit <- .dca_stan_surv_weibull(
    sample_size = nrow(surv_data),
    n_thr = n_thresholds,
    n_strategies = n_strategies,
    event_times_stacked = event_times_stacked,
    censored_times_stacked = censored_times_stacked,
    n_event_times_stacked = length(event_times_stacked),
    n_censored_times_stacked = length(censored_times_stacked),
    event_times_start_positions = event_times_start_positions,
    censored_times_start_positions = censored_times_start_positions,
    event_times_sizes = event_times_sizes,
    censored_times_sizes = censored_times_sizes,
    prediction_time = prediction_time,
    thresholds = thresholds,
    pos_post1 = posterior_positivity_pars$.shape1,
    pos_post2 = posterior_positivity_pars$.shape2,
    total_event_times = total_event_times,
    total_censored_times = total_censored_times,
    event_times_marginal = event_times_marginal,
    censored_times_marginal = censored_times_marginal,
    other_models_indices = other_models_indices,
    prior_only = as.numeric(prior_only),
    shape_prior = shape_prior,
    scale_prior = scale_prior,
    shape_prior_pars = shape_prior_pars,
    scale_prior_pars = scale_prior_pars,
    iter = iter,
    refresh = refresh,
    ...
  )

  dca_summary <- .extract_dca_surv_summary(
    fit = fit,
    summary_probs = summary_probs,
    thresholds = thresholds,
    strategies = strategies
  )

  output_data <- list(
    summary = dca_summary,
    thresholds = thresholds,
    .time = surv_data$.time,
    .status = surv_data$.status,
    .time_original = surv_data$.time * prediction_time,
    .data = .data,
    strategies = strategies,
    prediction_time = prediction_time,
    posterior_positivity_pars = posterior_positivity_pars
  )

  if (isTRUE(keep_fit)) {
    output_data[["fit"]] <- fit
  }

  if (isTRUE(keep_draws)) {
    output_data[["draws"]] <- .extract_dca_surv_draws(
      fit = fit,
      strategies = strategies
    )
  }

  .output <- structure(output_data, class = "BayesDCASurv")
  return(.output)
}

#' Extract posterior summaries for `BayesDCASurv` object
#' 
#' @param summary_probs Probabilities for posterior credible intervals (defaults to a 95% Cr.I.).
#' @param thresholds Decision thresholds -- within interval (0, 1).
#' @param strategies Vector of names of models or binary tests under assessment.
#' @return An object of class `BayesDCASurv`
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect str_c str_extract str_remove str_remove_all
#' @keywords internal
.extract_dca_surv_summary <- function(fit, # nolint
                                      summary_probs,
                                      thresholds,
                                      strategies) {
  .pars <- .get_relevant_pars()

  fit_summary <- rstan::summary(
    fit,
    pars = .pars,
    probs = summary_probs,
  )$summary %>%
    tibble::as_tibble(rownames = "par_name") %>%
    dplyr::select(-c(se_mean, sd, n_eff, Rhat)) %>% # nolint
    dplyr::rename(estimate := mean) %>% # nolint
    dplyr::mutate(
      # sorry about this terrible code but that's how life is sometimes
      threshold_ix = dplyr::case_when(
        str_detect(par_name, str_c(.pars[1:6], collapse = "|")) ~ str_extract(par_name, "\\[\\d+") %>% # nolint
          str_remove(string = ., pattern = "\\[") %>%
          as.integer(),
        str_detect(par_name, "positivity") ~ str_extract(par_name, "\\d+\\]") %>% # nolint
          str_remove(string = ., pattern = "\\]") %>%
          as.integer(),
        str_detect(par_name, "treat_all") ~ str_extract(par_name, "\\d+") %>% # nolint
          as.integer(),
        TRUE ~ NA_integer_
      ),
      decision_strategy_ix = dplyr::case_when(
        str_detect(par_name, str_c(.pars[1:6], collapse = "|")) ~ str_extract(par_name, "\\d+\\]") %>% # nolint
          str_remove(string = ., pattern = "\\]") %>%
          as.integer(),
        str_detect(par_name, str_c("positivity", collapse = "|")) ~ str_extract(par_name, "\\[\\d+") %>% # nolint
          str_remove(string = ., pattern = "\\[") %>%
          as.integer(),
        TRUE ~ NA_integer_
      ),
      threshold = ifelse(
        is.na(threshold_ix), # nolint
        NA_real_,
        thresholds[threshold_ix]
      ),
      decision_strategy_name = ifelse(
        is.na(decision_strategy_ix), # nolint
        NA_character_,
        strategies[decision_strategy_ix]
      ),
      par_name = stringr::str_extract(par_name, "\\w+")
    ) %>%
    dplyr::select(
      par_name, threshold, decision_strategy_name, # nolint
      dplyr::everything(), -dplyr::contains("ix")
    )

  # overall survival at prediction time
  overall_surv <- fit_summary %>%
    dplyr::filter(
      str_detect(par_name, "St_marginal") # nolint
    ) %>%
    dplyr::mutate(par_name = "overall_surv") %>%
    dplyr::select(-c(threshold, decision_strategy_name)) # nolint

  # positivity
  positivity <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "positivity")) # nolint

  # net benefit
  net_benefit <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "net_benefit")) # nolint

  # treat all (net benefit for treat all strategy)
  treat_all <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "treat_all")) %>% # nolint
    dplyr::select(-decision_strategy_name) # nolint

  # delta useful
  delta_useful <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "deta_useful")) # nolint

  # delta best
  delta_best <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "deta_best")) # nolint

  .summary <- structure(
    list(
      net_benefit = net_benefit,
      treat_all = treat_all,
      overall_surv = overall_surv,
      positivity = positivity,
      delta_useful = delta_useful,
      delta_best = delta_best
    ),
    class = "BayesDCASummary"
  )

  return(.summary)
}

#' Get relevant parameters to parse from Stan output (BayesDCASurv)
#' @keywords internal
.get_relevant_pars <- function() {
  return(
    c(
      "net_benefit",
      "delta_best", "delta_useful",
      "p_best", "p_useful",
      "best_competitor_nb",
      "positivity", "treat_all", "St_marginal"
    )
  )
}

#' @title Get posterior draws from time-to-event DCA stanfit
#'
#' @param fit A stanfit object.
#' @param strategies Vector of names of models
#' or binary tests under assessment.
#' @keywords internal
.extract_dca_surv_draws <- function(fit,
                                    strategies) {
  .pars <- .get_relevant_pars()

  stan_draws <- rstan::extract(
    fit,
    pars = .pars
  )

  .init_list <- function() {
    n_strategies <- length(strategies)
    x <- vector(mode = "list", length = n_strategies)
    setNames(x, strategies)
  }

  .draws <- list(
    overall_surv = as.vector(stan_draws$St_marginal),
    treat_all = stan_draws$treat_all,
    net_benefit = .init_list(),
    delta_useful = .init_list(),
    delta_best = .init_list(),
    p_useful = .init_list(),
    p_best = .init_list(),
    best_competitor_nb = .init_list()
  )

  for (i in seq_along(strategies)) {
    .name <- strategies[i]
    .draws[["net_benefit"]][[.name]] <- stan_draws$net_benefit[, , i]
    .draws[["delta_useful"]][[.name]] <- stan_draws$delta_useful[, , i]
    .draws[["delta_best"]][[.name]] <- stan_draws$delta_best[, , i]
    .draws[["p_useful"]][[.name]] <- stan_draws$p_useful[, , i]
    .draws[["p_best"]][[.name]] <- stan_draws$p_best[, , i]
    .draws[["best_competitor_nb"]][[.name]] <- stan_draws$best_competitor_nb[, , i] # nolint
  }

  return(.draws)
}


#' @title Print BayesDCASurv
#'
#' @param obj BayesDCASurv object
#' @export
print.BayesDCASurv <- function(obj, ...) {
  sample_size_info <- paste0(
    "N=", length(obj$.time), "\tevents=", sum(obj$.status)
  )
  os <- round(obj$summary$overall_surv[, -1] * 100, 1)
  cat(
    paste0(
      c(
        "BayesDCASurv\n",
        paste0("Number of thresholds: ", length(obj$thresholds)),
        sample_size_info,
        paste0("Prediction time point: ", obj$prediction_time),
        paste0(
          "Survival at prediction time: ",
          os$estimate, "% [",
          os$`2.5%`, "% \u2012 ",
          os$`97.5%`, "%]"
        ),
        "Models or tests: ",
        paste0(obj$strategies, collapse = ", ")
      ),
      collapse = "\n"
    )
  )
}

#' @title Get delta plot data for survival bayesDCA
#'
#' @param obj BayesDCA object
#' @param min_diff Minimum difference to call superior.
#' @param strategies Character vector with subset of decision strategies.
#' @param type One of "best", "useful", or "pairwise".
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @keywords internal
#' @return A ggplot object.
get_superiority_prob_plot_data_surv <- function(obj, # nolint: object_length_linter.
                                                min_diff = 0,
                                                strategies = NULL,
                                                type = c("best", "useful", "pairwise"),
                                                labels = NULL) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two strategies to plot pairwise comparison" = length(strategies) == 2 # nolint
    )
  }

  strategies <- validate_strategies(
    obj = obj, strategies = strategies
  )

  if (type == "pairwise") {
    nb1 <- obj$draws$net_benefit[[strategies[1]]]
    nb2 <- obj$draws$net_benefit[[strategies[2]]]
    .prob <- colMeans(nb1 - nb2 > min_diff)
    df <- tibble::tibble(
      prob = .prob,
      threshold = obj$thresholds
    )
  } else {
    df <- lapply(
      seq_along(strategies),
      function(j) {
        .m <- strategies[j]
        if (type == "best") {
          .prob <- colMeans(obj$draws$p_best[[.m]])
        } else {
          .prob <- colMeans(obj$draws$p_useful[[.m]])
        }
        tibble::tibble(
          prob = .prob,
          threshold = obj$thresholds,
          decision_strategy = .m,
          label = labels[.m]
        )
      }
    ) %>%
      dplyr::bind_rows()
  }
  return(df)
}

#' @title Get delta plot data for binary bayesDCA
#'
#' @param obj BayesDCA object
#' @param strategies Character vector with subset of decision strategies.
#' @param type One of "best", "useful", or "pairwise".
#' @param labels Named vector with label for each model or test.
#' @importFrom magrittr %>%
#' @keywords internal
#' @return A ggplot object.
get_delta_plot_data_surv <- function(obj,
                                     strategies = NULL,
                                     type = c("best", "useful", "pairwise"),
                                     labels = NULL) {
  type <- match.arg(type)
  if (type == "pairwise") {
    stopifnot(
      "Must specify two strategies to plot pairwise comparison" = length(strategies) == 2 # nolint
    )
  }

  strategies <- validate_strategies(
    obj = obj, strategies = strategies
  )

  if (is.null(obj$draws)) {
    msg <- "Retrieving posterior draws."
    message(cli::col_br_cyan(msg))
    if (!inherits(obj$fit, "stanfit")) {
      msg <- "FATAL - missing posterior draws and stanfit object. Recopute dca using keep_draws = TRUE."
    }
    obj$draws <- .extract_dca_surv_draws(
      fit = obj$fit,
      strategies = strategies
    )
  }

  if (type == "pairwise") {
    nb1 <- obj$draws$net_benefit[[strategies[1]]]
    nb2 <- obj$draws$net_benefit[[strategies[2]]]
    .delta <- nb1 - nb2
    qs <- matrixStats::colQuantiles(.delta, probs = c(.025, .975))
    df <- tibble::tibble(
      estimate = colMeans(.delta),
      `2.5%` = qs[, "2.5%"],
      `97.5%` = qs[, "97.5%"],
      threshold = obj$thresholds,
    )
  } else {
    df <- lapply(
      seq_along(strategies),
      function(i) {
        .m <- strategies[i]
        if (type == "best") {
          .delta <- obj$draws$delta_best[[.m]]
        } else {
          .delta <- obj$draws$delta_useful[[.m]]
        }
        qs <- matrixStats::colQuantiles(.delta, probs = c(.025, .975)) # nolint
        tibble::tibble(
          estimate = colMeans(.delta),
          `2.5%` = qs[, "2.5%"],
          `97.5%` = qs[, "97.5%"],
          threshold = obj$thresholds,
          decision_strategy = .m,
          label = labels[.m]
        )
      }
    ) %>%
      dplyr::bind_rows()
  }
  return(df)
}
