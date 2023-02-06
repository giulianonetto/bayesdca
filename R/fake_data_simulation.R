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
#' @keywords internal
#' @examples
#' d <- simulate_diagnostic_test_data(B = 2,
#'                                    N = 1000,
#'                                    true_p = 0.2,
#'                                    true_se = 0.9,
#'                                    true_sp = 0.9)
#' head(d)
simulate_diagnostic_test_data <- function(B = 100,  #nolint
                                          N = 500,  #nolint
                                          true_p = 0.2,
                                          true_se = 0.9,
                                          true_sp = 0.9,
                                          keep_true_pars = FALSE) {

  d <- rbinom(n = B, size = N, prob = true_p)
  tp <- rbinom(n = B, size = d, prob = true_se)
  tn <- rbinom(n = B, size = N - d, prob = true_sp)

  df <- data.frame(N, d, tp, tn)

  if (isTRUE(keep_true_pars)) {
    df <- cbind(df, true_p, true_se, true_sp)
  }

  return(df)
}


#' @title Simulate Prognostic Test Data
#' @importFrom magrittr %>%
#' @keywords internal
simulate_prognostic_model_data <- function(N = 200,
                                           .seed = 123) {
  stopifnot(require(simstudy))

  def <- defData(varname = "x1", formula = 0.5, dist = "binary")
  def <- defData(def, varname = "grp", formula = 0.5, dist = "binary")

  # Survival data definitions
  set.seed(.seed)
  sdef <- defSurv(varname = "survTime",
                  formula = "1.5*x1",
                  scale = "grp*50 + (1-grp)*25",
                  shape = "grp*1 + (1-grp)*1.5")
  sdef <- defSurv(sdef, varname = "censorTime", scale = 80, shape = 1)

  # Baseline data definitions
  .dtSurv <- genData(N, def) %>%  #nolint
    genSurv(
      sdef, timeName = "obsTime", censorName = "censorTime",
      eventName = "status", keepEvents = TRUE
    ) %>%
    dplyr::mutate(
      survTime = ifelse(survTime > 0, survTime, 0.001),  #nolint
      censorTime = ifelse(censorTime > 0, censorTime, 0.001),  #nolint
      obsTime = ifelse(status > 0, survTime, censorTime)  #nolint
    )

  return(.dtSurv)
}
