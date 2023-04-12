data{
  int<lower=0> T;
  real y0;
  vector[T] y;
}
parameters{
  real mu;                                // mean log volatility
  real < lower = 0, upper = 1 > phiT;     // persistence of volatility
  real < lower = 0 > s2;
  vector [T] h_std;                       // std log volatility time t
  vector < lower = 0 > [T]  l;
  real b0;
  real < lower = 0, upper = 1 > b1T;
  real b2;
  real <lower = 0 > v;
}
transformed parameters{
  real < lower = -1, upper = 1 > phi;
  real < lower = -1, upper = 1 > b1;
  real < lower = 0 > sigma;
  vector [T] h;                         // now h ~ normal(0, sigma)
  vector [T] mu_t;
  vector [T] sd_t;
  
  phi = 2 * phiT - 1;
  b1 = 2 * b1T - 1;
  sigma = sqrt(s2);
  
  //--- Volatilitys:
  h = h_std * sigma;
  h[1] /= sqrt(1 - phi * phi);         // rescale h[1]
  h += mu;
  for (t in 2:T) {
    h[t] += phi * (h[t - 1] - mu);
  }
  
  //--- Means:
  mu_t[1] = b0 + y0 * b1 + exp( h[1] ) * b2;
  for(t in 2:T){
    mu_t[t] = b0 + y[t-1] * b1 + exp( h[t] ) * b2;
  }
  //--- Variances:
  for(t in 1:T){
    sd_t[t] = exp( h[t] / 2 ) / sqrt( l[t] );
  }
  
}
model{
  // --- prioris log - volatiliti
  mu ~ normal(0, 10);
  phiT ~ beta(20, 1.5);
  s2 ~ inv_gamma(2.5, 0.025);
  
  //--- prioris component means:
  b0 ~ normal(0, 10);
  b1T ~ beta(5, 1.5);
  b2 ~ normal(0, 10);
  
  //--- prioris component heavy-tailed:
  v ~ gamma(2, 2);
  for(t in 1:T){
    l[t] ~ gamma(0.5 * v, 0.5 * v);  
  }
  
  /*
  //--- Sampling volatilitys:
  //h_std ~ std_normal();
  h[1] ~ normal(mu, sigma / sqrt(1 - phi * phi));
  for (t in 2:T){
    h[t] ~ normal(mu + phi * ( h[t - 1] -  mu ), sigma);
  }
  */
  
  //--- Sampling observations:
  y ~ normal( mu_t, sd_t );
  /*
  y[1] ~ normal(b0 + y0 * b1 + exp( h[1] ) * b2,  exp( h[1]/2 ) );
  for(t in 2:T){
    y[t] ~ normal(b0 + y0 * b1 + exp( h[t] ) * b2,  exp( h[t]/2 ) );
  }
  */
}
