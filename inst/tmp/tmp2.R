library(tidyverse)
library(rstan)

thresholds <- seq(0, 0.9, .1)
names(thresholds) <- thresholds
df_pop <- read_tsv("inst/tmp/example.tsv")
set.seed(5851)
# ix <- sample(seq_len(nrow(df_pop)), 137)
ix <- sample(seq_len(nrow(df_pop)), 1e3)
df_pop <- df_pop[ix, ] %>%
  mutate(
    outcomes = survival::Surv(obsTime, status)
  ) %>%
  select(outcomes, p_hat)
surv_data <- data.frame(
  .time = unname(df_pop[["outcomes"]][, 1]), # observed time-to-event
  .status = unname(df_pop[["outcomes"]][, 2]) # 1 if event, 0 if censored
)
posterior_positivity_pars <- bayesDCA:::get_positivity_posterior_parameters(
  .prediction_data = df_pop[, -1],
  .thresholds = thresholds
)
# cutpoints <- c(0, 0.4, 0.6, 0.8, 1)
cutpoints <- bayesDCA:::get_cutpoints(1, surv_data$.time[surv_data$.status == 1])
posterior_surv_pars0 <- bayesDCA:::get_survival_posterior_parameters(
  .prediction_data = NA,
  .surv_data = surv_data,
  .models_or_tests = 1,
  .cutpoints = cutpoints,
  .thresholds = 0,
  .prior_scaling_factor = 1 / 50,
  .prior_means = NULL
)
n_models_or_tests <- 1
deaths_and_times <- get_death_pseudo_counts(
  .prediction_data = df_pop[, -1],
  .surv_data = surv_data,
  .cutpoints = cutpoints,
  .models_or_tests = names(df_pop[, -1]),
  .thresholds = thresholds
)
.data <- list(
  n_thr = length(thresholds),
  n_models = 1,
  n_intervals = length(cutpoints),
  thresholds = thresholds,
  time_exposed = deaths_and_times$total_exposure_times, # bayesDCA:::get_survival_time_exposed(1, cutpoints),
  time_exposed_prediction = get_survival_time_exposed2(1, cutpoints),
  deaths = deaths_and_times$death_pseudo_counts, # n_models-long list with n_intervals rows and n_thr columns
  pos_post1 = posterior_positivity_pars$.shape1,
  pos_post2 = posterior_positivity_pars$.shape2,
  posterior_alpha0 = posterior_surv_pars0$.alpha,
  posterior_beta0 = posterior_surv_pars0$.beta,
  other_models_indices = lapply(
    1:n_models_or_tests,
    function(i) seq_len(n_models_or_tests)[-i]
  ),
  prior_only = 0
)
model <- rstan::stan_model("inst/tmp/dca_surv_hierarchical.stan")
fit <- sampling(model, data = .data, refresh = 1, cores = 4, iter = 4000)
s <- rstan::summary(fit)$summary %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("par_name")

d <- s %>%
  as_tibble() %>%
  filter(str_detect(par_name, "net_benefit")) %>%
  mutate(
    i = str_extract(par_name, "\\d+") %>% as.numeric(),
    thr = thresholds[i]
  )
ggplot(d, aes(thr)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3) +
  geom_line(aes(y = mean)) +
  geom_line(
    data = true_nb %>% filter(.thr %in% d$thr),
    aes(.thr, nb),
    inherit.aes = FALSE,
    linetype = "longdash",
    color = "red"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(-1, 1))

left_join(d, true_nb %>%
  rename(thr = .thr), by = "thr") %>%
  ggplot(aes(mean, nb)) +
  geom_point() +
  geom_abline()
