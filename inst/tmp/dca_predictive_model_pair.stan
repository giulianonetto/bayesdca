data {
  
  // Input Data
  int<lower = 1> n_thr;  // number of thresholds
  vector<lower = 0>[n_thr] N;  // sample sizes for model 1
  vector<lower = 0>[n_thr] tp_m1; // true positives for model 1
  vector<lower = 0>[n_thr] tp_m2; // true positives for model 2
  vector<lower = 0>[n_thr] tn_m1; // true negatives for model 1
  vector<lower = 0>[n_thr] tn_m2; // true negatives for model 2
  vector<lower = 0>[n_thr] d;  // diseased persons or cases for model 1
  vector<lower = 0, upper = 1>[n_thr] thresholds;
  // Prior parameters for prevalence
  real<lower = 0> prior_p1;
  real<lower = 0> prior_p2;
  // Prior parameters for Sensitivity
  vector<lower = 0>[n_thr] prior_Se1_m1;
  vector<lower = 0>[n_thr] prior_Se2_m1;
  vector<lower = 0>[n_thr] prior_Se1_m2;
  vector<lower = 0>[n_thr] prior_Se2_m2;
  // Prior parameters for Specificity
  vector<lower = 0>[n_thr] prior_Sp1_m1;
  vector<lower = 0>[n_thr] prior_Sp2_m1;
  vector<lower = 0>[n_thr] prior_Sp1_m2;
  vector<lower = 0>[n_thr] prior_Sp2_m2;
}
transformed data {
  // Posterior parameters for prevalence
  real<lower = 0> post_p1;
  real<lower = 0> post_p2;
  // Posterior parameters for Sensitivity
  vector<lower = 0>[n_thr] post_Se1_m1;
  vector<lower = 0>[n_thr] post_Se2_m1;
  vector<lower = 0>[n_thr] post_Se1_m2;
  vector<lower = 0>[n_thr] post_Se2_m2;
  // Posterior parameters for Specificity
  vector<lower = 0>[n_thr] post_Sp1_m1;
  vector<lower = 0>[n_thr] post_Sp2_m1;
  vector<lower = 0>[n_thr] post_Sp1_m2;
  vector<lower = 0>[n_thr] post_Sp2_m2;
  // prevalence
  post_p1 = d[1] + prior_p1;
  post_p2 = N[1] - d[1] + prior_p2;
  // sensitivity for model 1
  post_Se1_m1 = tp_m1 + prior_Se1_m1;
  post_Se2_m1 = d - tp_m1 + prior_Se2_m1;
  // sensitivity for model 2
  post_Se1_m2 = tp_m2 + prior_Se1_m2;
  post_Se2_m2 = d - tp_m2 + prior_Se2_m2;
  // specificity for model 1
  post_Sp1_m1 = tn_m1 + prior_Sp1_m1;
  post_Sp2_m1 = N - d - tn_m1 + prior_Sp2_m1;
  // specificity for model 2
  post_Sp1_m2 = tn_m2 + prior_Sp1_m2;
  post_Sp2_m2 = N - d - tn_m2 + prior_Sp2_m2;
}
parameters {
  real<lower=0, upper = 1> p;
  vector<lower=0, upper = 1>[n_thr] Se_m1;
  vector<lower=0, upper = 1>[n_thr] Se_m2;
  vector<lower=0, upper = 1>[n_thr] Sp_m1;
  vector<lower=0, upper = 1>[n_thr] Sp_m2;
}

transformed parameters {
  vector<lower = 0, upper = 1>[n_thr] p_vec;
  for (i in 1:n_thr) {
    p_vec[i] = p;
  }
}
model {
  // Analytical posterior for prevalence (or outcome proportion)
  p ~ beta(post_p1, post_p2);
  // Analytical posteriors for model 1
  Se_m1 ~ beta(post_Se1_m1, post_Se2_m1);
  Sp_m1 ~ beta(post_Sp1_m1, post_Sp2_m1);
  // Analytical posteriors for model 2
  Se_m2 ~ beta(post_Se1_m2, post_Se2_m2);
  Sp_m2 ~ beta(post_Sp1_m2, post_Sp2_m2);
}
generated quantities {
  
  // Treat all (same for both models)
  vector[n_thr] treat_all;
  // Net benefir calculation for model 1
  vector[n_thr] net_benefit_m1;
  vector[n_thr] delta_m1;
  // Net benefir calculation for model 2
  vector[n_thr] net_benefit_m2;
  vector[n_thr] delta_m2;
  // Pairwise comparison between model 1 and model 2
  vector[n_thr] delta_m1_vs_m2;
  vector[n_thr] prob_better_m1_vs_m2;
  for (i in 1:n_thr) {
    // treat all calculation for threhsold i
    treat_all[i] = 1*p_vec[i] - (1-p_vec[i])*(1-0)*(thresholds[i]/(1-thresholds[i]));
    // single threshold calculation for model 1
    net_benefit_m1[i] = Se_m1[i]*p_vec[i] - (1-p_vec[i])*(1-Sp_m1[i])*(thresholds[i]/(1-thresholds[i]));
    if (treat_all[i] > 0) {
      delta_m1[i] = net_benefit_m1[i] - treat_all[i];
    } else {
      delta_m1[i] = net_benefit_m1[i];
    }
    // single threshold calculation for model 2
    net_benefit_m2[i] = Se_m2[i]*p_vec[i] - (1-p_vec[i])*(1-Sp_m2[i])*(thresholds[i]/(1-thresholds[i]));
    if (treat_all[i] > 0) {
      delta_m2[i] = net_benefit_m2[i] - treat_all[i];
    } else {
      delta_m2[i] = net_benefit_m2[i];
    }
    // Pairwise comparison: model 1 versus model 2
    delta_m1_vs_m2[i] = net_benefit_m1[i] - net_benefit_m2[i];
    prob_better_m1_vs_m2[i] = delta_m1_vs_m2[i] > 0;
  }
}
