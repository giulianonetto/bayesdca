#' Predictive Model Evaluation Data.
#'
#' A dataset containing the outcomes (0 or 1), probability predictions
#' for an example model, and results from an example binary test.
#' Binary test generated using sensitivity
#' and specificity of 80%. Based on `rmda::dcaData`.
#'
#' @format A data frame with 500 rows and 3 variables:
#' \describe{
#'   \item{outcomes}{binary outcomes (0 or 1, int.)}
#'   \item{predictions}{predicted probabilities (num.)}
#'   \item{binary_test}{binary test results (0 or 1, int.)}
#'   ...
#' }
#' @source \url{https://giulianonetto.github.io/bayesdca/}
"PredModelData"

#' Survival Data for Bayesian DCA.
#'
#' A dataset containing the the time-to-event `outcomes` (a `survival::Surv` object),
#' the predictions from a prognostic model, and the results from a binary
#' prognostic test. Based on `dcurves::df_surv`. The intended time horizon for prediction
#' is exactly 1.
#'
#' @format A data frame with 500 rows and 3 variables:
#' \describe{
#'   \item{outcomes}{survival outcomes (`survival::Surv` object)}
#'   \item{model_predictions}{predicted probabilities of event (at time 1).}
#'   \item{binary_test}{binary prognostic test results (0 or 1, int.)}
#'   ...
#' }
#' @source \url{https://giulianonetto.github.io/bayesdca/}
"dca_survival_data"

#' Survival data from Docking, 2021 - Figures 3B-3E
#'
#' A dataset containing the survival data from [Docking, 2021](https://www.nature.com/articles/s41467-021-22625-y).
#' Dataset includes selected data used for Figures 3B to 3E.
#' Observations with missing survival time (in weeks)
#' or missing survival events were excluded.
#' AML prognostic score (APS) values were kept as an example of continuous biomarker.
#'
#' @format A data frame with 635 rows and 5 variables:
#' \describe{
#'   \item{fig_id}{Character indicating corresponding figure in the original paper.}
#'   \item{aps_lasso_score}{APS value}
#'   \item{aps_lasso_class}{APS tercile computed within study cohort.}
#'   \item{survival_event}{Death is 1, censored/slive is 0.}
#'   \item{survival_weeks}{Overall survival in weeks.}
#'   ...
#' }
#' @source \url{https://giulianonetto.github.io/bayesdca/}
"Docking2021Data"
