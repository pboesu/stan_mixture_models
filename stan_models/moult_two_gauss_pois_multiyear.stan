/* two gaussian curves with a poisson error model */
data {
  int<lower=0> N; //number of observations
  int<lower=0> Y; //number of seasons
  int<lower=0> counts[Y,N]; 
  real times[Y,N]; 
  real min_times;
  //real max_times;
} 
parameters {
 
  real total_birds[Y]; 
  real<lower=min_times,upper=175> mean1[Y];
  real<lower=0> diff_mean2[Y];
  real<lower=0,upper=0.5> proportion[Y];
  real<lower=0> sigma1[Y]; 
  real<lower=0> sigma2[Y];
  //real<lower=0> sigma_total;
  //real intercept; // intercept of regression model
  //real beta_date; // regression coefficient for return date effect
  //real<lower=0> lm_sigma; //linear model sigma
} 
transformed parameters {
  real m[Y,N];
  //real post_total[Y];
  //real proportion = 0.5;
  for(j in 1:Y){
  for (i in 1:N){ 
    m[j,i] = total_birds[j] * proportion[j] / sqrt(2*pi()*sigma1[j]^2) * exp( -0.5 * (times[j,i] - mean1[j])^2 / sigma1[j]^2) + total_birds[j]*(1-proportion[j]) / sqrt(2*pi()*sigma2[j]^2) * exp( -0.5 * (times[j,i] - (mean1[j] + diff_mean2[j]))^2 / sigma2[j]^2);
  }
  //post_total[j] = (1-proportion[j])*total_birds[j];
  }
} 
model {
  //vector[Y] pred_counts;
  
  for (j in 1:Y){
  counts[j,] ~ poisson(m[j,]); 
  sum(counts[j,]) ~ poisson(total_birds[j]); //likelihood on the total number of birds. might need a neg bin here?
    //regression model
  //pred_counts[j] = intercept + beta_date * mean1[j];
  //post_total[j] ~ normal(pred_counts[j], lm_sigma);
  }

  
  
  
  total_birds ~ gamma(.001, .001);
  //sigma_total ~ cauchy(0,1);
  proportion ~ beta(2, 25); 
  mean1 ~ normal(162, 2); //moult arrrival dates from williams & croxall 1991 interpreted as 95% quantiles
  diff_mean2 ~ normal(24, 1); //moult duration from williams & croxall 1991
  sigma1 ~ cauchy(0, 0.1); 
  sigma2 ~ cauchy(0, 0.1); 
  //intercept ~ gamma(0.001,0.001);
  //beta_date ~ normal(0,100);
  //lm_sigma ~ cauchy(0,10);
}
generated quantities{
  real mean2[Y];
  real post_total[Y];
  //vector[N] log_lik;
  for (j in 1:Y){
  mean2[j] = mean1[j] + diff_mean2[j];
  post_total[j] = (1-proportion[j])*total_birds[j];
  }
  //for (n in 1:N) log_lik[n] = poisson_lpmf(counts[n]| m[n]); 
}
