data{
  int<lower=0> T;
  real y0;
  vector[T] y;
}
parameters{
  real mu;                         // mean log volatility
  real<lower=-1,upper=1> phiT;     // persistence of volatility
  real<lower=0> s2;
  vector[T] h_std;  // std log volatility time t
  real<lower=-1,upper=1> b1T;
  real b0;
  real b2;
}
transformed parameters{
  real<lower=-1,upper=1> phi;
  real<lower=-1,upper=1> b1;
  real<lower=0> sigma;
  vector[T] h;  // now h ~ normal(0, sigma)
  vector[T] mu_t;
  phi = (2*phiT - 1);
  b1 = (2*b1T - 1);
  sigma = sqrt(s2);
  
  //--- Volatilitys:
  h = h_std * sigma;
  h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
  h += mu;
  for (t in 2:T) {
    h[t] += phi * (h[t - 1] - mu);
  }
  
  //--- Means:
  mu_t[1] = b0 + y0*b1 + exp( h[1] ) * b2;
  for(t in 2:T){
    mu_t[t] = b0 + b1 * y[t-1]  + exp( h[t] ) * b2;
  }
}
model{
  // --- prioris
  mu ~ normal(0, 3.2);
  phiT ~ beta(20, 1.5);
  b0 ~ normal(0, 3.2);
  b1T ~ beta(5, 1.5);
  b2 ~ normal(0, 3.2);
  s2 ~ inv_gamma(2.5, 0.025);
  
  //--- Sampling volatilitys:
  h_std ~ std_normal();
  
  //--- Sampling observations:
  y ~ normal( mu_t, exp(h/2) );
}
