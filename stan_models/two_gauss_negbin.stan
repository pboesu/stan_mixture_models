data {
int<lower=0> N; 
int<lower=0> counts[N]; 
real times[N]; 
} 
parameters {
real total_birds; 
real<lower=54> mean1;
real diff_mean2;
real<lower=0,upper=1> proportion;
real<lower=0> sigma1; 
real<lower=0> sigma2;
real<lower=0> phi;
} 
transformed parameters {
real m[N];
for (i in 1:N) 
m[i] = total_birds * proportion / sqrt(2*pi()*sigma1^2) * exp( -0.5 * (times[i] - mean1)^2 / sigma1^2) + total_birds*(1-proportion) / sqrt(2*pi()*sigma2^2) * exp( -0.5 * (times[i] - (mean1 + diff_mean2))^2 / sigma2^2);
} 
model {

counts ~ neg_binomial_2(m, phi); 

total_birds ~ gamma(.001, .001); 
proportion ~ beta(200, 200); 
mean1 ~ gamma(2, 2); 
diff_mean2 ~ gamma(2, 2);
sigma1 ~ cauchy(0, 1); 
sigma2 ~ cauchy(0, 1); 
phi ~ cauchy(0,1);
}
generated quantities{
real mean2;
vector[N] log_lik;
mean2 = mean1 + diff_mean2;
for (n in 1:N) log_lik[n] = neg_binomial_2_lpmf(counts[n]| m[n], phi); 
}
