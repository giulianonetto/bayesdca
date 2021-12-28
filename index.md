
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
68%. Decision Curve Analysis helps us find an answer.

# Installation

You can install the development version of bayesDCA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("giulianonetto/bayesdca")
```

# Examples

## Clinical prediction model

``` r
library(bayesDCA)
data(PredModelData)
fit <- dca_predictive_model(df = PredModelData)
plot(fit)
```

## Diagnostic test (binary)

``` r
library(bayesDCA)
fit <- dca_diagnostic_test(N = 500, d = 83, tp = 77, tn = 378)
plot(fit)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Diagnostic test (continuous)

Under construction.

## Compare Net Benefit

``` r
library(bayesDCA)
fit1 <- dca_diagnostic_test(N = 500, d = 100, tp = 70, tn = 380)
fit2 <- dca_diagnostic_test(N = 500, d = 100, tp = 90, tn = 200)
compare_net_benefit("Test A" = fit1, "Test B" = fit2)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
