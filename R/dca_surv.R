#' Fit Bayesian Decision Curve Analysis using Stan for survival outcomes
#'
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#' @keywords internal
.dca_stan_surv <- function(n_thr,
                           n_models_or_tests,
                           n_intervals,
                           thresholds,
                           time_exposed,
                           posterior_alpha,
                           posterior_beta,
                           posterior_alpha0,
                           posterior_beta0,
                           pos_post1,
                           pos_post2,
                           other_models_indices,
                           iter = 4000,
                           refresh = 0, ...) {
  thresholds <- pmin(thresholds, 0.9999) # odds(1) = Inf

  standata <- list(
    n_thr = n_thr,
    n_models = n_models_or_tests,
    n_intervals = n_intervals,
    thresholds = thresholds,
    time_exposed = time_exposed,
    posterior_alpha = posterior_alpha,
    posterior_beta = posterior_beta,
    posterior_alpha0 = posterior_alpha0,
    posterior_beta0 = posterior_beta0,
    pos_post1 = pos_post1,
    pos_post2 = pos_post2,
    refresh = refresh
  )

  dots <- list(...)
  if ("control" %in% names(dots)) {
    control <- dots[["control"]]
  } else {
    control <- list(adapt_delta = 0.9)
  }

  .model <- stanmodels$dca_time_to_event
  stanfit <- rstan::sampling(.model,
    data = standata,
    control = control,
    iter = iter,
    refresh = refresh, ...
  )
  return(stanfit)
}

#' Fit Bayesian Decision Curve Analysis using Stan for survival outcomes (hierarchical model)
#'
#' @param refresh Control verbosity of [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
#' @param ... Arguments passed to [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains).
#' @return An object of class [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.html) returned by [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html) (e.g. iter, chains)
#' @keywords internal
.dca_stan_surv2 <- function(n_thr,
                            n_models_or_tests,
                            n_intervals,
                            thresholds,
                            time_exposed_prediction,
                            deaths,
                            time_exposed,
                            posterior_alpha0,
                            posterior_beta0,
                            prior_only,
                            pos_post1,
                            pos_post2,
                            other_models_indices,
                            iter = 4000,
                            refresh = 0, ...) {
  thresholds <- pmin(thresholds, 0.9999) # nolint odds(1) = Inf

  standata <- list(
    n_thr = n_thr,
    n_models = n_models_or_tests,
    n_intervals = n_intervals,
    thresholds = thresholds,
    time_exposed = time_exposed,
    time_exposed_prediction = time_exposed_prediction,
    deaths = deaths,
    pos_post1 = pos_post1,
    pos_post2 = pos_post2,
    posterior_alpha0 = posterior_alpha0,
    posterior_beta0 = posterior_beta0,
    other_models_indices = other_models_indices,
    prior_only = prior_only,
    refresh = refresh
  )

  dots <- list(...)
  if ("control" %in% names(dots)) {
    control <- dots[["control"]]
  } else {
    control <- list(adapt_delta = 0.9)
  }

  .model <- stanmodels$dca_time_to_event2
  stanfit <- rstan::sampling(
    .model,
    data = standata,
    control = control,
    iter = iter,
    refresh = refresh, ...
  )
  return(stanfit)
}

#' Bayesian Decision Curve Analysis for Predictive Models and Binary tests
#'
#' @export
#' @return An object of class `BayesDCASurv`
#' @importFrom magrittr %>%
#' @examples
#' data(dca_survival_data)
#' fit <- dca_surv(dca_survival_data, prediction_time = 1, cores = 4)
#' plot(fit)
dca_surv <- function(.data,
                     prediction_time,
                     thresholds = seq(0, 0.5, 0.02),
                     keep_draws = TRUE,
                     keep_fit = FALSE,
                     summary_probs = c(0.025, 0.975),
                     cutpoints = NULL,
                     prior_scaling_factor = 1 / 3,
                     prior_means = NULL,
                     iter = 4000,
                     refresh = 0,
                     ...) {
  if (colnames(.data)[1] != "outcomes") {
    stop("Missing 'outcomes' column as the first column in input .data")
  }

  stopifnot(
    "'outcomes' column must be a Surv object. " = survival::is.Surv(.data[["outcomes"]])
  )

  # avoid thresholds in {0, 1}
  thresholds <- thresholds %>%
    pmin(0.99) %>%
    pmax(1e-9)

  # preprocess .data
  model_or_test_names <- colnames(.data)[-1]
  prediction_data <- data.frame(.data[, -1])
  colnames(prediction_data) <- model_or_test_names
  surv_data <- data.frame(
    .time = unname(.data[["outcomes"]][, 1]), # observed time-to-event
    .status = unname(.data[["outcomes"]][, 2]) # 1 if event, 0 if censored
  )

  event_times <- surv_data$.time[surv_data$.status > 0]

  # preprocess cutpoints
  if (is.null(cutpoints)) {
    cutpoints <- get_cutpoints(prediction_time, event_times)
  }

  # make sure zero is included and cutpoints are correctly ordered
  cutpoints <- sort(unique(c(0, cutpoints)))
  # stopifnot(
  #   "prediction_time must be at most max(cutpoints)" = prediction_time <= max(cutpoints)
  # )


  events_per_interval <- get_events_per_interval(cutpoints, event_times)
  epi_text <- paste(
    names(events_per_interval), events_per_interval,
    sep = ": ", collapse = "  "
  ) %>%
    stringr::str_replace("Inf]", "Inf)")
  if (!(all(events_per_interval >= min_events_per_interval()))) {
    msg <- paste0(
      "Please updade cutpoints, too few events per interval (at least ",
      min_events_per_interval(), " required):\n",
      epi_text
    )
    stop(msg)
  }

  msg <- paste0(
    "Survival estimation with the following intervals (total event counts):\n",
    epi_text,
    collapse = "\n"
  )
  message(msg)

  n_models_or_tests <- ncol(prediction_data)
  n_thresholds <- length(thresholds)
  n_intervals <- length(cutpoints)
  time_exposed <- get_survival_time_exposed(prediction_time, cutpoints)

  posterior_surv_pars <- get_survival_posterior_parameters(
    .prediction_data = prediction_data,
    .surv_data = surv_data,
    .models_or_tests = colnames(prediction_data),
    .cutpoints = cutpoints,
    .thresholds = thresholds,
    .prior_scaling_factor = prior_scaling_factor,
    .prior_means = prior_means
  )
  posterior_surv_pars0 <- get_survival_posterior_parameters(
    .prediction_data = NA,
    .surv_data = surv_data,
    .models_or_tests = 1,
    .cutpoints = cutpoints,
    .thresholds = 0,
    .prior_scaling_factor = prior_scaling_factor,
    .prior_means = prior_means
  )
  posterior_positivity_pars <- get_positivity_posterior_parameters(
    .prediction_data = prediction_data,
    .thresholds = thresholds
  )

  other_models_indices <- lapply(
    1:n_models_or_tests,
    function(i) (1:n_models_or_tests)[-i]
  )

  fit <- .dca_stan_surv(
    n_thr = n_thresholds,
    n_models = n_models_or_tests,
    n_intervals = n_intervals,
    thresholds = thresholds,
    time_exposed = time_exposed,
    posterior_alpha = posterior_surv_pars$.alpha,
    posterior_beta = posterior_surv_pars$.beta,
    posterior_alpha0 = posterior_surv_pars0$.alpha,
    posterior_beta0 = posterior_surv_pars0$.beta,
    pos_post1 = posterior_positivity_pars$.shape1,
    pos_post2 = posterior_positivity_pars$.shape2,
    other_models_indices = other_models_indices,
    iter = iter,
    refresh = refresh,
    ...
  )

  dca_summary <- .extract_dca_surv_summary(
    fit = fit,
    summary_probs = summary_probs,
    thresholds = thresholds,
    model_or_test_names = model_or_test_names
  )

  output_data <- list(
    summary = dca_summary,
    thresholds = thresholds,
    .time = surv_data$.time,
    .status = surv_data$.status,
    .data = .data,
    cutpoints = cutpoints,
    model_or_test_names = model_or_test_names,
    prediction_time = prediction_time,
    posterior_surv_pars = posterior_surv_pars,
    posterior_positivity_pars = posterior_positivity_pars
  )

  if (isTRUE(keep_fit)) {
    output_data[["fit"]] <- fit
  }

  if (isTRUE(keep_draws)) {
    output_data[["draws"]] <- .extract_dca_surv_draws(
      fit = fit,
      model_or_test_names = model_or_test_names
    )
  }

  .output <- structure(output_data, class = "BayesDCASurv")
  return(.output)
}


#' Bayesian Decision Curve Analysis for Survival outcomes (hierarchical)
#'
#' @export
#' @return An object of class `BayesDCASurv`
#' @importFrom magrittr %>%
#' @examples
#' data(dca_survival_data)
#' fit <- dca_surv(dca_survival_data, prediction_time = 1, cores = 4)
#' plot(fit)
dca_surv2 <- function(.data,
                      prediction_time,
                      thresholds = seq(0, 0.5, 0.02),
                      keep_draws = TRUE,
                      keep_fit = FALSE,
                      summary_probs = c(0.025, 0.975),
                      cutpoints = NULL,
                      prior_scaling_factor = 1 / 50,
                      prior_means = NULL,
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

  # avoid thresholds in {0, 1}
  thresholds <- thresholds %>%
    pmin(0.99) %>%
    pmax(1e-9)

  # preprocess .data
  model_or_test_names <- colnames(.data)[-1]
  prediction_data <- data.frame(.data[, -1])
  colnames(prediction_data) <- model_or_test_names
  surv_data <- data.frame(
    .time = unname(.data[["outcomes"]][, 1]), # observed time-to-event
    .status = unname(.data[["outcomes"]][, 2]) # 1 if event, 0 if censored
  )

  event_times <- surv_data$.time[surv_data$.status > 0]

  # preprocess cutpoints
  if (is.null(cutpoints)) {
    cutpoints <- get_cutpoints(prediction_time, event_times)
  }

  # make sure zero is included and cutpoints are correctly ordered
  cutpoints <- sort(unique(c(0, cutpoints)))

  events_per_interval <- get_events_per_interval(cutpoints, event_times)
  epi_text <- paste(
    names(events_per_interval), events_per_interval,
    sep = ": ", collapse = "  "
  ) %>%
    stringr::str_replace("Inf]", "Inf)")
  if (!(all(events_per_interval >= min_events_per_interval()))) {
    msg <- paste0(
      "Please updade cutpoints, too few events per interval (at least ",
      min_events_per_interval(), " required):\n",
      epi_text
    )
    stop(msg)
  }

  msg <- paste0(
    "Survival estimation with the following intervals (total event counts):\n",
    epi_text,
    collapse = "\n"
  )
  message(msg)

  n_models_or_tests <- ncol(prediction_data)
  n_thresholds <- length(thresholds)
  n_intervals <- length(cutpoints)
  time_exposed_prediction <- get_survival_time_exposed(
    .prediction_time = prediction_time,
    .cutpoints = cutpoints
  )

  posterior_surv_pars0 <- get_survival_posterior_parameters(
    .prediction_data = NA,
    .surv_data = surv_data,
    .models_or_tests = 1,
    .cutpoints = cutpoints,
    .thresholds = 0,
    .prior_scaling_factor = prior_scaling_factor,
    .prior_means = prior_means
  )
  posterior_positivity_pars <- get_positivity_posterior_parameters(
    .prediction_data = prediction_data,
    .thresholds = thresholds
  )

  other_models_indices <- lapply(
    1:n_models_or_tests,
    function(i) (1:n_models_or_tests)[-i]
  )

  deaths_and_times <- get_death_pseudo_counts(
    .prediction_data = prediction_data,
    .surv_data = surv_data,
    .cutpoints = cutpoints,
    .models_or_tests = model_or_test_names,
    .thresholds = thresholds
  )

  fit <- .dca_stan_surv2(
    n_thr = n_thresholds,
    n_models = n_models_or_tests,
    n_intervals = n_intervals,
    thresholds = thresholds,
    time_exposed_prediction = time_exposed_prediction,
    time_exposed = deaths_and_times$total_exposure_times,
    deaths = deaths_and_times$death_pseudo_counts,
    posterior_alpha0 = posterior_surv_pars0$.alpha,
    posterior_beta0 = posterior_surv_pars0$.beta,
    pos_post1 = posterior_positivity_pars$.shape1,
    pos_post2 = posterior_positivity_pars$.shape2,
    other_models_indices = other_models_indices,
    prior_only = as.numeric(prior_only),
    iter = iter,
    refresh = refresh,
    ...
  )

  dca_summary <- .extract_dca_surv_summary(
    fit = fit,
    summary_probs = summary_probs,
    thresholds = thresholds,
    model_or_test_names = model_or_test_names
  )

  output_data <- list(
    summary = dca_summary,
    thresholds = thresholds,
    .time = surv_data$.time,
    .status = surv_data$.status,
    .data = .data,
    deaths_and_times = deaths_and_times,
    cutpoints = cutpoints,
    model_or_test_names = model_or_test_names,
    prediction_time = prediction_time,
    posterior_surv_pars0 = posterior_surv_pars0,
    posterior_positivity_pars = posterior_positivity_pars
  )

  if (isTRUE(keep_fit)) {
    output_data[["fit"]] <- fit
  }

  if (isTRUE(keep_draws)) {
    output_data[["draws"]] <- .extract_dca_surv_draws(
      fit = fit,
      model_or_test_names = model_or_test_names
    )
  }

  .output <- structure(output_data, class = "BayesDCASurv")
  return(.output)
}

#' Extract posterior summaries for `BayesDCASurv` object
#'
#' @return An object of class `BayesDCASurv`
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect str_c str_extract str_remove str_remove_all
#' @keywords internal
.extract_dca_surv_summary <- function(fit, # nolint
                                      summary_probs,
                                      thresholds,
                                      model_or_test_names) {
  .pars <- c(
    "net_benefit", "delta", "prob_better_than_soc",
    "highest_nb_other_than_model_j",
    "positivity", "treat_all", "St_marginal"
  )

  fit_summary <- rstan::summary(
    fit,
    pars = .pars,
    probs = summary_probs,
  )$summary %>%
    tibble::as_tibble(rownames = "par_name") %>%
    dplyr::select(-c(se_mean, sd, n_eff, Rhat)) %>% # nolint
    dplyr::rename(estimate := mean) %>% # nolint
    dplyr::mutate(
      threshold_ix = dplyr::case_when(
        str_detect(par_name, str_c(.pars[1:4], collapse = "|")) ~ str_extract(par_name, "\\[\\d+") %>% # nolint
          str_remove(string = ., pattern = "\\[") %>%
          as.integer(),
        str_detect(par_name, str_c(.pars[5:6], collapse = "|")) ~ str_extract(par_name, "\\d+\\]") %>% # nolint
          str_remove(string = ., pattern = "\\]") %>%
          as.integer(),
        TRUE ~ NA_integer_
      ),
      model_or_test_ix = dplyr::case_when(
        str_detect(par_name, str_c(.pars[1:4], collapse = "|")) ~ str_extract(par_name, "\\d+\\]") %>% # nolint
          str_remove(string = ., pattern = "\\]") %>%
          as.integer(),
        str_detect(par_name, str_c(.pars[5], collapse = "|")) ~ str_extract(par_name, "\\[\\d+") %>% # nolint
          str_remove(string = ., pattern = "\\[") %>%
          as.integer(),
        TRUE ~ NA_integer_
      ),
      threshold = ifelse(
        is.na(threshold_ix), # nolint
        NA_real_,
        thresholds[threshold_ix]
      ),
      model_or_test_name = ifelse(
        is.na(model_or_test_ix), # nolint
        NA_character_,
        model_or_test_names[model_or_test_ix]
      ),
      par_name = stringr::str_extract(par_name, "\\w+"),
      par_name = ifelse(
        str_detect(
          par_name,
          "highest_nb_other_than_model_j"
        ),
        "best_competitor_nb",
        par_name
      )
    ) %>%
    dplyr::select(
      par_name, threshold, model_or_test_name, # nolint
      dplyr::everything(), -dplyr::contains("ix")
    )

  # overall survival at prediction time
  overall_surv <- fit_summary %>%
    dplyr::filter(
      str_detect(par_name, "St_marginal") # nolint
    ) %>%
    dplyr::mutate(par_name = "overall_surv") %>%
    dplyr::select(-c(threshold, model_or_test_name)) # nolint

  # positivity
  positivity <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "positivity")) # nolint

  # net benefit
  net_benefit <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "net_benefit")) # nolint

  # treat all (net benefit for treat all strategy)
  treat_all <- fit_summary %>%
    dplyr::filter(str_detect(par_name, "treat_all")) %>% # nolint
    dplyr::select(-model_or_test_name) # nolint

  .summary <- structure(
    list(
      net_benefit = net_benefit,
      treat_all = treat_all,
      overall_surv = overall_surv,
      positivity = positivity
    ),
    class = "BayesDCASummary"
  )

  return(.summary)
}

#' @title Get posterior draws from time-to-event DCA stanfit
#'
#' @param fit A stanfit object.
#' @param model_or_test_names Vector of names of models
#' or binary tests under assessment.
#' @keywords internal
.extract_dca_surv_draws <- function(fit,
                                    model_or_test_names) {
  .pars <- c(
    "net_benefit", "delta",
    "prob_better_than_soc",
    "highest_nb_other_than_model_j",
    "positivity", "treat_all", "St_marginal"
  )

  stan_draws <- rstan::extract(
    fit,
    pars = .pars
  )

  .draws <- list(
    overall_surv = stan_draws$St_marginal %>% as.vector(),
    treat_all = stan_draws$treat_all,
    net_benefit = list(),
    delta_default = list(),
    prob_best = list()
  )

  for (i in seq_along(model_or_test_names)) {
    .name <- model_or_test_names[i]
    .draws[["net_benefit"]][[.name]] <- stan_draws$net_benefit[, , i]
    .draws[["delta_default"]][[.name]] <- stan_draws$delta[, , i]
    .draws[["prob_best"]][[.name]] <- stan_draws$prob_better_than_soc[, , i]
    .draws[["best_competitor_nb"]][[.name]] <- stan_draws$highest_nb_other_than_model_j[, , i] # nolint
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
        paste0(obj$model_or_test_names, collapse = ", ")
      ),
      collapse = "\n"
    )
  )
}
