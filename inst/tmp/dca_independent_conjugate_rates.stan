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
  int<lower = 0> fp;
  real<lower = 0> post_p1;
  real<lower = 0> post_p2;
  real<lower = 0> post_tpr1;
  real<lower = 0> post_tpr2;
  real<lower = 0> post_fpr1;
  real<lower = 0> post_fpr2;
  // Observed false-positives
  fp = N - d - tn;
  // TPR
  post_tpr1 = tp + prior_Se1;
  post_tpr2 = N - tp + prior_Se2;
  // FPR
  post_fpr1 = fp + prior_Sp1;
  post_fpr2 = N - d - fp + prior_Sp1;
  // prevalence
  post_p1 = d + prior_p1;
  post_p2 = N - d + prior_p2;
}
parameters {
  real<lower=0, upper = 1> p;
  real<lower=0, upper = 1> tpr;
  real<lower=0, upper = 1> fpr;
}
model {
  // Analytical posterior
  p ~ beta(post_p1, post_p2);
  tpr ~ beta(post_tpr1, post_tpr2);
  fpr ~ beta(post_fpr1, post_fpr2);
}
generated quantities{

  vector[n_thresholds] net_benefit;
  vector[n_thresholds] treat_all;
  for (i in 1:n_thresholds) {
    net_benefit[i] = tpr - fpr*(thresholds[i]/(1-thresholds[i]));
    treat_all[i] = 1*p - (1-p)*(1-0)*(thresholds[i]/(1-thresholds[i]));
  }

}

