data {

  // Input Data
  int<lower = 1> B;  // number of data sets (N, d, tp, tn, thresholds)
  vector<lower = 0>[B] N;
  vector<lower = 0>[B] tp;
  vector<lower = 0>[B] tn;
  vector<lower = 0>[B] d;
  vector<lower = 0, upper = 1>[B] thresholds;
  // Priors
  real<lower = 0> prior_p1;
  real<lower = 0> prior_p2;
  vector<lower = 0>[B] prior_Se1;
  vector<lower = 0>[B] prior_Se2;
  vector<lower = 0>[B] prior_Sp1;
  vector<lower = 0>[B] prior_Sp2;
}
transformed data {
  real<lower = 0> post_p1;
  real<lower = 0> post_p2;
  vector<lower = 0>[B] post_Se1;
  vector<lower = 0>[B] post_Se2;
  vector<lower = 0>[B] post_Sp1;
  vector<lower = 0>[B] post_Sp2;
  // prevalence
  post_p1 = d[1] + prior_p1;
  post_p2 = N[1] - d[1] + prior_p2;
  // sensitivity
  post_Se1 = tp + prior_Se1;
  post_Se2 = d - tp + prior_Se2;
  // specificity
  post_Sp1 = tn + prior_Sp1;
  post_Sp2 = N - d - tn + prior_Sp2;
}
parameters {
  real<lower=0, upper = 1> p;
  vector<lower=0, upper = 1>[B] Se;
  vector<lower=0, upper = 1>[B] Sp;
}

transformed parameters {
  vector<lower = 0, upper = 1>[B] p_vec;
  for (i in 1:B) {
    p_vec[i] = p;
  }
}
model {
  // Analytical posterior
  p ~ beta(post_p1, post_p2);
  Se ~ beta(post_Se1, post_Se2);
  Sp ~ beta(post_Sp1, post_Sp2);
}
generated quantities{

  vector[B] net_benefit;
  vector[B] treat_all;
  for (i in 1:B) {
    net_benefit[i] = Se[i]*p_vec[i] - (1-p_vec[i])*(1-Sp[i])*(thresholds[i]/(1-thresholds[i]));
    treat_all[i] = 1*p_vec[i] - (1-p_vec[i])*(1-0)*(thresholds[i]/(1-thresholds[i]));
  }

}
