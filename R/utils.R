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
