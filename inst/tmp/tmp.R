
priors <- expand.grid(
  alphas = c(1, 2, 2.5, 5, 10),
  sigmas = c(1, 2, 2.5, 5, 10)
)

surv_plots <- vector("list", nrow(priors))
par(mfrow = c(5, 5))
for (i in seq_len(nrow(priors))) {
  a <- priors[i, "alphas"]
  s <- priors[i, "sigmas"]
  cat(paste0("alpha sd ", a, " sigma sd ", s))
  .standata <- list(
    n_events = sum(d$status),
    n_censored = sum(d$status == 0L),
    event_times = d$obsTime[d$status == 1L],
    censored_times = d$obsTime[d$status == 0L],
    pred_time = pred_time,
    prior_scale = s,
    prior_sd = a,
    prior_only = 1
  )
  .stanfit <- sampling(
    m2,
    data = .standata,
    chains = 2, iter = 5000,
    cores = 2, refresh = 0
  )
  .draws <- tidybayes::tidy_draws(.stanfit) %>%
    select(.draw, alpha, sigma) %>%
    mutate(
      posterior_surv = purrr::pmap(
        list(alpha, sigma),
        function(a, s) {
          tibble(ts = times, s = exp(-(times / s)^a))
        }
      )
    )

  .lab <- paste0(
    "alpha=", a, " sigma=", s
  )
  hist(.draws$sigma, breaks = 200, xlim = c(0, 50), main = .lab)
  hist(.draws$alpha, breaks = 200, xlim = c(0, 50), add = TRUE, col = "red")
  surv_plots[[i]] <- .draws %>%
    unnest(posterior_surv) %>%
    filter(near(round(ts, 3), 1, .006)) %>%
    ggplot(aes(alpha, sigma, color = s)) +
    geom_point() +
    coord_cartesian(ylim = c(0, 20), xlim = c(0, 20)) +
    facet_wrap(~ts) +
    scale_color_viridis_b(breaks = seq(0, 1, .2)) +
    theme(legend.key.height = unit(2, "cm")) +
    geom_rug()
}



system.time({
  x <- dca_surv_weibull(
    d,
    prediction_time = 12,
    iter = 3000,
    cores = 2,
    chains = 2,
    refresh = 1,
    thresholds = seq(0, .9, length = 51),
    mean_mu = 0,
    sd_mu = 5,
    mean_log_alpha = 1,
    sd_log_alpha = 1,
    prior_only = FALSE,
    keep_fit = TRUE
  )
})
