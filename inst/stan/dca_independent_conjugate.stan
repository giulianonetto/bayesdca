data {
  // Input Data
  int<lower = 0> N;
  int<lower = 0> tp;
  int<lower = 0> tn;
  int<lower = 0> d;
  int<lower = 0> n_thresholds;
  vector[n_thresholds] thresholds;
  // Priors
  real<lower = 0> prior_p1;
  real<lower = 0> prior_p2;
  real<lower = 0> prior_Se1;
  real<lower = 0> prior_Se2;
  real<lower = 0> prior_Sp1;
  real<lower = 0> prior_Sp2;
}
transformed data {
  real<lower = 0> post_p1;
  real<lower = 0> post_p2;
  real<lower = 0> post_Se1;
  real<lower = 0> post_Se2;
  real<lower = 0> post_Sp1;
  real<lower = 0> post_Sp2;
  // prevalence
  post_p1 = d + prior_p1;
  post_p2 = N - d + prior_p2;
  // sensitivity
  post_Se1 = tp + prior_Se1;
  post_Se2 = d - tp + prior_Se2;
  // specificity
  post_Sp1 = tn + prior_Sp1;
  post_Sp2 = N - d - tn + prior_Sp2;
}
parameters {
  real<lower=0, upper = 1> p;
  real<lower=0, upper = 1> Se;
  real<lower=0, upper = 1> Sp;
}
model {
  // Analytical posterior
  p ~ beta(post_p1, post_p2);
  Se ~ beta(post_Se1, post_Se2);
  Sp ~ beta(post_Sp1, post_Sp2);
}
generated quantities{
  vector[n_thresholds] net_benefit;
  vector[n_thresholds] treat_all;
  vector[n_thresholds] delta;
  for (i in 1:n_thresholds) {
    net_benefit[i] = Se*p - (1-p)*(1-Sp)*(thresholds[i]/(1-thresholds[i]));
    treat_all[i] = 1*p - (1-p)*(1-0)*(thresholds[i]/(1-thresholds[i]));
    if (treat_all[i] > 0) {
      delta[i] = net_benefit[i] - treat_all[i];
    } else {
      delta[i] = net_benefit[i];
    }
  }
}
