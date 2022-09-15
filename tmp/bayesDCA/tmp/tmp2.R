library(rstan)
library(tidyverse)

thresholds = seq(0, .5, .01)
names(thresholds) <- thresholds
model <- rstan::stan_model("inst/stan/dca.stan")
.data <- list(
  N = 1e3, d = 100, tp = 95, tn = 810,
  n_thresholds = length(thresholds), thresholds = thresholds,
  prior_p1 = 1, prior_p2 = 1/100,
  prior_Se1 = 1, prior_Se2 = 1,
  prior_Sp1 = 1, prior_Sp2 = 1
)
fit <- sampling(model, data=.data, refresh=0)
s <- rstan::summary(fit)$summary %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("par_name")

true_nb <- map(.data$thresholds, ~ {
  tibble(nb = .9*.1 - .1*.9*(.x/(1-.x)))
}) %>% bind_rows(.id = "thr") %>%
  mutate(thr = as.numeric(thr))

d <- s %>% as_tibble() %>%
  filter(str_detect(par_name, "net_benefit")) %>%
  mutate(
    i = str_extract(par_name, "\\d+") %>% as.numeric(),
    thr = thresholds[i]
  )
ggplot(d, aes(thr)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3) +
  geom_line(aes(y=mean)) +
  geom_line(
    data = true_nb,
    aes(thr, nb),
    inherit.aes = FALSE,
    linetype = "longdash",
    color = "red"
  ) +
  geom_hline(yintercept = c(0,.1),
             linetype = "dashed") +
  theme_bw()

# DCA pair fitting
model <- rstan::stan_model("inst/tmp/dca_predictive_model_pair.stan")
b <- sample(c(-1, 1), nrow(PredModelData), replace = T)
PredModelData2 <- PredModelData %>% dplyr::mutate(
  predictions = plogis(log(predictions/(1-predictions)) + rnorm(nrow(.), sd=.5))
  )
thresholds = seq(0, .5, .01)
names(thresholds) <- thresholds
m1_data <- get_thr_data(PredModelData$outcomes,
                        PredModelData$predictions,
                        thresholds = thresholds)
m2_data <- get_thr_data(PredModelData2$outcomes,
                        PredModelData2$predictions,
                        thresholds = thresholds)
n_thr <- length(thresholds)
.data <- list(
  n_thr = n_thr,
  N = m1_data$N,
  tp_m1 = m1_data$tp,
  tp_m2 = m2_data$tp,
  tn_m1 = m1_data$tn,
  tn_m2 = m2_data$tn,
  d = m1_data$d,
  thresholds = thresholds,
  prior_p1 = 1,
  prior_p2 = 1,
  prior_Se1_m1 = rep(1, n_thr),
  prior_Se2_m1 = rep(1, n_thr),
  prior_Se1_m2 = rep(1, n_thr),
  prior_Se2_m2 = rep(1, n_thr),
  prior_Sp1_m1 = rep(1, n_thr),
  prior_Sp2_m1 = rep(1, n_thr),
  prior_Sp1_m2 = rep(1, n_thr),
  prior_Sp2_m2 = rep(1, n_thr)
)
fit <- sampling(model, data=.data, refresh=0, cores = 4)
s <- rstan::summary(fit)$summary %>%
  data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("par_name")
d <- rstan::extract(fit)
j = 20; plot(d$net_benefit_m1[,j], d$net_benefit_m2[,j], main = thresholds[j]);  cor(d$net_benefit_m1[,j], d$net_benefit_m2[,j])
plot(
  thresholds,
  sapply(1:51, function(i) {
    cor(
      d$net_benefit_m1[,i],
      d$net_benefit_m2[,i]
    )
  })
)

# DCA list fitting

model2 <- rstan::stan_model("inst/tmp/dca_predictive_model_list.stan")

.data2 <- list(
  n_thr = n_thr,
  n_models = 2,
  N = m1_data$N,
  d = m1_data$d,
  tp = cbind(m1_data$tp, m2_data$tp),
  tn = cbind(m1_data$tn, m2_data$tn),
  thresholds = thresholds,
  prior_p1 = 1,
  prior_p2 = 1,
  prior_Se1 = cbind(rep(1, n_thr), rep(1, n_thr)),
  prior_Se2 = cbind(rep(1, n_thr), rep(1, n_thr)),
  prior_Sp1 = cbind(rep(1, n_thr), rep(1, n_thr)),
  prior_Sp2 = cbind(rep(1, n_thr), rep(1, n_thr))
)
fit2 <- sampling(model2, data=.data2, cores = 4)
s <- rstan::summary(fit)$summary %>%
  data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("par_name")
d <- rstan::extract(fit)
j = 20; plot(d$net_benefit_m1[,j], d$net_benefit_m2[,j], main = thresholds[j]);  cor(d$net_benefit_m1[,j], d$net_benefit_m2[,j])
plot(
  thresholds,
  sapply(1:51, function(i) {
    cor(
      d$net_benefit_m1[,i],
      d$net_benefit_m2[,i]
    )
  })
)
