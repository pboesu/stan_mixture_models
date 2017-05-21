data {
 int<lower = 0> N;
 vector[N] y;
 real min_y;
}

parameters {
  real<lower=min_y> mu1;
  real<lower=0> diff_mu;
  real<lower=0> sigma[2];
  real<lower=0, upper=1> theta;
}

model {
 sigma ~ normal(0, 2);
 mu1 ~ normal(60, 2);
 diff_mu ~ normal(12, 0.5);//from Williams 1991 
 theta ~ beta(10, 10);
 for (n in 1:N)
   target += log_mix(theta,
                     normal_lpdf(y[n] | mu1, sigma[1]),
                     normal_lpdf(y[n] | mu1+diff_mu, sigma[2]));
}
