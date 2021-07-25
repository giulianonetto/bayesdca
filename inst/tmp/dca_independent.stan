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
parameters {
  real<lower=0, upper = 1> p;
  real<lower=0, upper = 1> Se;
  real<lower=0, upper = 1> Sp;
}
model {
  // Priors
  p ~ beta(prior_p1, prior_p2);
  Se ~ beta(prior_Se1, prior_Se2);
  Sp ~ beta(prior_Sp1, prior_Sp2);

  // Model
  d ~ binomial(N, p);
  tp ~ binomial(d, Se);
  tn ~ binomial(N-d, Sp);
}
generated quantities{

  vector[n_thresholds] net_benefit;
  vector[n_thresholds] treat_all;
  for (i in 1:n_thresholds) {
    net_benefit[i] = Se*p - (1-p)*(1-Sp)*(thresholds[i]/(1-thresholds[i]));
    treat_all[i] = 1*p - (1-p)*(1-0)*(thresholds[i]/(1-thresholds[i]));
  }

}

