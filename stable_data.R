############################################################################
###################### alpha - stable data
library( stabledist )

mu = -1
phi = 0.95
sigma = 0.13
b0 = 0.001
b1 = 0.11
b2 = -0.21
y0 = 0
# er
a = 1.92    # a e ( 0, 2 ] 
T = 3e3
y = h = l = rep(0, T)

data_seed = sample( 1:1e6, 1 )
set.seed( data_seed )
for(t in 1:T){
  if(t == 1){
    h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
    y[t] = rstable(n = 1, 
                   alpha = a, 
                   beta = 0, 
                   gamma = exp(h[t]/2),
                   delta = b0 + b1 * y0 + b2 * exp( h[t] ) 
                   )
  }else{
    h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
    y[t] = rstable(n = 1, 
                   alpha = a, 
                   beta = 0, 
                   gamma = exp(h[t]/2),
                   delta = b0 + b1 * y[t-1] + b2 * exp( h[t] ) 
    )
  }
}

par( mfrow = c( 1, 2) )
plot( y, type = 'l', main = mean(y) )
hist( y, breaks = 40 )
par( mfrow = c(1, 1) )


