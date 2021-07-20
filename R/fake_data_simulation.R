#' @title Simulate Diagnostic Test Data
#'
#' @param B Integer indicating number of simulated datasets.
#' @param N Integer indicating sample size in each dataset.
#' @param true_p Proportion indicating true prevalence in each dataset.
#' @param true_se Proportion indicating true sensitivity in each dataset.
#' @param true_sp Proportion indicating true specificity in each dataset.
#' @param keep_true_pars Logical indicating whether to keep true
#' parameter values in the output dataframe.
#' @return A Bx7 data frame with simulated data and provided parameters.
#' @details Diagnostic test data is simulated according to the following
#' mechanism: for each of the `B` simulations, `d` *diseased* persons
#' are sampled from a binomial of size `N` and parameter `true_p`. Then,
#' `tp` *true positive* outcomes are sampled from a binomial of size `d`
#' and parameter `true_se`, and `tn` *true negative* outcomes are
#' sampled from a binomial of size `N-d` with parameter `true_sp`. This
#' simulates `B` diagnostic accuracy studies in which a number of
#' diseased and non-diseased persons are sampled from a population;
#' the expected proportion of diseased persons correctly detected as
#' such is the sensitivity, and the number of non-diseased persons
#' correctly detected as such is the specificity.
#' @importFrom magrittr %>%
#' @export
simulate_diagnostic_test_data <- function(B = 100,
                                          N = 500,
                                          true_p = 0.2,
                                          true_se = 0.9,
                                          true_sp = 0.9,
                                          keep_true_pars = TRUE) {

  d <- rbinom(n = B, size = N, prob = true_p)
  tp <- rbinom(n = B, size = d, prob = true_se)
  tn <- rbinom(n = B, size = N - d, prob = true_sp)

  df <- data.frame(
    N, d, tp, tn
  )

  if (isTRUE(keep_pars)) {
    df <- cbind(df, true_p, true_se, true_sp)
  }

  return(df)
}
