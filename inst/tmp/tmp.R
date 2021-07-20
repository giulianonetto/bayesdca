library(rstan)
library(tidyverse)

thresholds = seq(0, .5, .01)
names(thresholds) <- thresholds
model <- stan_model("inst/stan/dca.stan")
.data <- list(
  N = 1e3, d = 100, tp = 95, tn = 810,
  # n_thresholds = length(thresholds), thresholds = thresholds
  n_thresholds = 1, thresholds = array(.2)
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

