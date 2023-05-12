############################################################################
################################# SVM-SMN Model ############################
# y_t = b0 + b1*y_{t-1} + b2*e(h_t) + e_t/lambda_t^{1/2}
# h_t = mu + phi(h_{t-1} - mu) + sigma*n_{t}
############################################################################
###################### Normal: l ~ 1 

mu = -1
phi = 0.985
sigma = 0.13
b0 = 0.08
b1 = 0.15
b2 = - 0.23
y0 = 0
T = 2e3
y = h = rep(0, T)

#seed = sample(1:1e6, 1)
set.seed( 788284 )

for(t in 1:T){
  if(t == 1){
    h[t] = rnorm(1, mean = mu, sd = sigma*sqrt(1/(1-phi**2)))
    y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp(h[t]), exp(h[t]/2) )
  }else{
    h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu) ), sd = sigma)
    y[t] = rnorm(1, b0 + b1*y[t - 1] + b2 * exp(h[t]), exp(h[t]/2) )
  }
}

par(mfrow = c(1, 2))
plot(y, type = 'l', main = 'log - retornos')
hist(y, breaks = 40, main = 'histograma')
par(mfrow = c(1, 1))


