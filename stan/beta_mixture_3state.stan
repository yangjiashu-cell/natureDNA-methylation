// 3-component Beta mixture for one sample's fCpG beta values (0–1).
// Aligns with Nature 2025 Methods (Stan Beta-mixture discretization) in structure;
// full equivalence requires the same priors / MCMC settings as the original analysis.
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  simplex[3] theta;
  vector<lower=0.05, upper=80>[3] alpha;
  vector<lower=0.05, upper=80>[3] beta;
}
model {
  theta ~ dirichlet(rep_vector(1, 3));
  alpha ~ lognormal(0, 1.5);
  beta ~ lognormal(0, 1.5);
  for (n in 1 : N) {
    vector[3] log_p = log(theta);
    for (k in 1 : 3) {
      log_p[k] = log_p[k] + beta_lpdf(y[n] | alpha[k], beta[k]);
    }
    target += log_sum_exp(log_p);
  }
}
