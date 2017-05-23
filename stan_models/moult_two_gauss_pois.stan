/* two gaussian curves with a poisson error model */
data {
  int<lower=0> N; 
  int<lower=0> counts[N]; 
  real times[N]; 
  real min_times;
} 
parameters {
 
  real total_birds; 
  real<lower=min_times> mean1;
  real<lower=0> diff_mean2;
  real<lower=0,upper=0.5> proportion;
  real<lower=0> sigma1; 
  real<lower=0> sigma2;
  //real<lower=0> sigma_total;
} 
transformed parameters {
  real m[N];
  //real proportion = 0.5;
  for (i in 1:N) 
    m[i] = total_birds * proportion / sqrt(2*pi()*sigma1^2) * exp( -0.5 * (times[i] - mean1)^2 / sigma1^2) + total_birds*(1-proportion) / sqrt(2*pi()*sigma2^2) * exp( -0.5 * (times[i] - (mean1 + diff_mean2))^2 / sigma2^2);
} 
model {
  
  
  counts ~ poisson(m); 
  sum(counts) ~ poisson(total_birds); //likelihood on the total number of birds. might need a neg bin here?
  
  total_birds ~ gamma(.001, .001);
  //sigma_total ~ cauchy(0,1);
  proportion ~ beta(2, 25); 
  mean1 ~ normal(162, 2); //moult arrrival dates from williams & croxall 1991 interpreted as 95% quantiles
  diff_mean2 ~ normal(24, 1); //moult duration from williams & croxall 1991
  sigma1 ~ cauchy(0, 0.1); 
  sigma2 ~ cauchy(0, 0.1); 
}
generated quantities{
  real mean2;
  vector[N] log_lik;
  mean2 = mean1 + diff_mean2;
  for (n in 1:N) log_lik[n] = poisson_lpmf(counts[n]| m[n]); 
}
