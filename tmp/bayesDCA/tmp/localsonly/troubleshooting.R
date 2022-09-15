q <- matrixStats::colQuantiles(fit1$draws$delta, probs = c(.025, .5, .975))
q <- tibble(
  estimate = q[, "50%"],
  lower = q[, "2.5%"],
  upper = q[, "97.5%"],
  thr = fit1$thresholds
)
q2 <- matrixStats::colQuantiles(fit1$draws$net_benefit - fit1$draws$treat_all,
                                probs = c(.025, .5, .975))

q2 <- tibble(
  estimate = q2[, "50%"],
  lower = q2[, "2.5%"],
  upper = q2[, "97.5%"],
  thr = fit1$thresholds
)

plot(df[df$x=='nb', 'y'], df[df$x=='ta', 'y'], xlab = "nb", ylab = 'ta')
abline(0,1,col='red')
data.frame(
  x = rep(c("nb", 'ta'), each=4e3),
  y = c(
    fit1$draws$net_benefit[,8] ,
    fit1$draws$treat_all[,8]
  )
) %>%
  ggplot(aes(y, fill = x)) +
  geom_density(alpha = .3) + theme_bw()
