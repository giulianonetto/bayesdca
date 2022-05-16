
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesDCA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Perform Bayesian Decision Curve Analysis for clinical prediction models
and diagnostic tests.

When validating a clinical prediction model, you may end up with AUC
0.79 and a slight miscalibration. How do you know if that is good enough
for your model to be clinically useful? The same question can be asked
if you have a diagnostic test with sensitivity of 75% and specificity of
68%, for instance. Decision Curve Analysis helps us find an answer - see
[Vickers, van Calster & Steyerberg,
2019](https://diagnprognres.biomedcentral.com/articles/10.1186/s41512-019-0064-7)
for an introduction to DCA. Here, we use Bayesian methods to accurately
quantify uncertainty in our decisions curves - powered by
[Stan](https://mc-stan.org/).

# Installation

You can install the development version of bayesDCA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("giulianonetto/bayesdca")
```

# Running Bayesian DCA

You can use `bayesDCA` to evaluate predictive models as well as binary
tests. See the [last section](#continuous) for information on continuous
tests.

All you need is a `data.frame` with a column named `outcomes` (0 or 1)
and one column for each model or test being evaluated. In the example
below, the `PredModelData` includes the probability predictions from a
model and the results from a binary test. The names of these columns
don’t matter (except for the `outcomes` column).

``` r
library(bayesDCA)
data(PredModelData)
head(PredModelData)
#>   outcomes predictions binary_test
#> 1        0  0.01280653           0
#> 2        0  0.13981948           0
#> 3        0  0.03566458           0
#> 4        0  0.02351731           0
#> 5        0  0.00863298           0
#> 6        0  0.00959754           0
```

We set `cores = 4` to speed up MCMC sampling with
[Stan](https://mc-stan.org/).

``` r
fit <- dca(PredModelData, cores = 4)
plot(fit)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Comparing two decision strategies

Say you want to infer whether the predictions from the model yield a
better decision strategy than the binary test – i.e., you want to
compare their decision curves. Then:

``` r
compare_dca(fit, 
            models_or_tests = c("predictions", "binary_test"))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Comparing against default stratefies

Default strategies include treating all patients and treating no
patients. If the decision threshold is above the prevalence of the
disease, treating no patient is better than treating all patients. If
you run `compare_dca` and specify only one decision strategy, `bayesDCA`
will compare this strategy against the the appropriate default for each
decision threshold.

``` r
compare_dca(fit, 
            models_or_tests = "predictions")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## Continuous tests, scores, gene signatures, etc.

[Categorization of continuous predictors is strongly
discouraged](https://www.prognosisresearch.com/videos-categorisation),
and continuous tests are no different. This also applies to prognostic
or diagnostic scores, expression signatures, and all of the like. As
decisions are often categorical, it’s better to categorize outcome
probabilities and then map probability thresholds back to corresponding
test values. You end up with a test cutoff either way, but only one is
clinically motivated - see [Myth
2](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-019-1425-3#Sec3)
in [Wynants et al,
2019](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-019-1425-3).
This risk stratification procedure can be applied seamlessly by fitting
single-predictor models (e.g. `Disease ~ test_value`). This leads to
clinically-informed test cutoffs as opposed to observed-data
arbitrariness (e.g. categorizing by quantile of the score, which yields
arbitrary risk distributions in the score strata).

Notice that choosing an appropriate risk threshold for your continuous
test properly balances your requirements for sensitivity and specificity
according to the clinical context - the cost of each correct or
incorrect decision. Whereas one may wish to maximize sensitivity and
specificity simultaneously (e.g. pick the “elbow” in a ROC curve), that
procedure effectively assumes that false positives have the same costs
as false negatives, which is usually dramatically incorrect. For
instance, a screening test surely should penalize false negatives way
more heavily than false positives, while for a confirmatory test the
cost of false positives would be more important. Hence, the resulting
test cutoff should be optimized by clinical reasoning, not by
data-driven artifacts.

Once you have your test cutoff, you can compute `tp` and `tn` and
proceed with `dca_binary_test` as usual. Of course, a likely better
approach would be to use your model to take advantage of the
(non-linear) relationship between test values and outcome probability,
but that turns your diagnostic/prognostic test study into a predictive
model study - you would use `dca_predictive_model` in that case.

Future versions of `bayesDCA` should include an easy API for this task.
