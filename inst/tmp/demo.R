library(bayesDCA)
library(tidyverse)
data("PredModelData")

PredModelData %>% glimpse()
fit1 <- dca_predictive_model(outcomes = PredModelData$outcomes,
                             predictions = PredModelData$predictions)
plot(fit1)

# same prevalence as in PredModelData (12%)
fit2 <- dca_binary_test(N = 500, d = 60, tp = 54, tn = 365)
fit3 <- dca_binary_test(N = 500, d = 60, tp = 58, tn = 400)
plot_dca_list("model" = fit1, "test A" = fit2, "test B" = fit3)
compare_dca(model = fit1, "test A" = fit2)
compare_dca(testB = fit3, model = fit1)
compare_dca(testB = fit3, testA = fit2)
