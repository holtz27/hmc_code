############################################################################
################################# SVM-SMN Model ############################
# y_t = b0 + b1 * y_{t-1} + b2 * e(h_t) + e_t/lambda_t^{1/2}
# h_t = mu + phi * ( h_{t-1} - mu ) + sigma * n_{t}
############################################################################
###################### t-Student: lambda|v ~ Gama(v/2, v/2) 

mu = -1.0
phi = 0.985
sigma = 0.13
b0 = 0.01
b1 = 0.1
b2 = -0.23
y0 = 0
v = 17.89 #( log(v) )
T = 2e3
y = h = l = rep(0, T)

#data_seed = sample( 1:1e6, 1 )
# 634323
set.seed( 634323 )
for(t in 1:T){
  if(t == 1){
    l[t] = rgamma(1, shape = v/2, rate = v/2)
    h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
    y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp( h[t] ), exp(h[t]/2) / sqrt( l[t] ))
  }else{
    l[t] = rgamma(1, shape = v/2, rate = v/2)
    h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
    y[t] = rnorm(1, b0 + b1 * y[t-1] + b2 * exp(h[t]), exp(h[t]/2) / sqrt( l[t] ))
  }
}

par( mfrow = c( 1, 2) )
plot( y, type = 'l', main = 'log-retornos' )
hist( y, breaks = 40 )
par( mfrow = c(1, 1) )

moments::kurtosis(y)


