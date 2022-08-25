
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
if you have, for instance, a binary diagnostic test with sensitivity of
75% and specificity of 68%. Decision Curve Analysis helps us find an
answer - see [Vickers, van Calster & Steyerberg,
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
tests.

All you need is a `data.frame` with a column named `outcomes` (0 or 1)
and one column for each model or test being evaluated. In the example
below, the `PredModelData` includes the probability predictions from a
model (`"predictions" column`) and the results from a binary test
(`"binary_test"` column). The names of these columns don’t matter
(except for the `outcomes` column, which should always be present).

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

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" height="50%" />

## Comparing two decision strategies

Say you want to infer whether the predictions from the model yield a
better decision strategy than the binary test – i.e., you want to
compare their decision curves. Then:

``` r
compare_dca(fit)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" height="50%" />

## Comparing against default strategies

Default strategies include treating all patients and treating no
patients. If you run `compare_dca` and specify only one decision
strategy, `bayesDCA` will compare this strategy against the the
appropriate default for each decision threshold.

``` r
compare_dca(fit, models_or_tests = "binary_test")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" height="50%" />

## Using external information to estimate prevalence

Say you are validating tests using a nested case-control study, so the
prevalence parameter must come from the larger sample from which cases
and controls were selected. Another example is when you want to use an
external estimate of prevalence (say from a large prospective
cross-sectional study). You can do so by passing the
`external_prevalence_data` argument to `dca`. Notice that you may want
to adjust the default `thresholds` argument as well: in the example
below, the prevalence is around 60/20,000 = 0.3%, so we use pretty low
thresholds.

``` r
fit <- dca(PredModelData, cores = 4,
           external_prevalence_data = c(60,20000),
           thresholds = seq(0, 0.01, 0.001))
plot(fit)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" height="50%" />

Notice also that the `external_prevalence_data` information is used only
to estimate the prevalence, which is then used in the net benefit
calculation. Sensitivity and specificity use the original information in
`PredModelData`.
