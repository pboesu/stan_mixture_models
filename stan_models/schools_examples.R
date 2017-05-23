schools_data <- list(
  J = 8,
  y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
)

library(rstan)

model1 <- "data {
  int<lower=0> J;          // number of schools 
  real y[J];               // estimated treatment effects
  real<lower=0> sigma[J];  // s.e. of effect estimates 
}
parameters {
  real mu; 
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta;
  theta = mu + tau * eta;
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}
"

fit1 <- stan(
  model_code = model1,  # Stan program
  data = schools_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4,              # number of cores (using 2 just for the vignette)
  refresh = 1000          # show progress every 'refresh' iterations
)


model2 <- "data {
int<lower=0> J; // # schools
real y[J]; // estimated treatment
real<lower=0> sigma[J]; // std err of effect
}
parameters {
real theta[J]; // school effect
real mu; // mean for schools
real<lower=0> tau; // variance between schools
}
model {
theta ~ normal(mu, tau);
y ~ normal(theta, sigma);
}
"

fit2 <- stan(
  model_code = model2,  # Stan program
  data = schools_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4           # number of cores (using 2 just for the vignette)
          # show progress every 'refresh' iterations
)

traceplot(fit2)
plot(fit2, pars=c("theta"))


plot(fit1, pars = "theta")

get_posterior_overall_mean <- function(stanfit, pars){
  get_posterior_mean(stanfit, pars)[,length(get_seeds(stanfit))+1]
}

plot(get_posterior_overall_mean(fit1, "theta"), get_posterior_overall_mean(fit2, "theta"))
abline(0,1,col="red")
