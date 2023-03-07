path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/Parâmetros/theta/jacobiano/theta.cpp'
Rcpp::sourceCpp(path)

#theta = (mu, w, gama)
#param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)

seed = sample(1:1e3, 1)

N = 1e3
z = rmhmc(N, eps = c(0.015, 0.015, 0.015), 
          min_L = 10, max_L = 15, 
          theta_init = c( -1, atanh( -0.95 ), log( 0.25 ) ), 
          param = list(y = c(y0, y), 
                       h = rnorm(T), 
                       l = rgamma(T, 2), 
                       mu_0 = 0, s_0 = 3.2, 
                       a_phi = 20, b_phi = 1.5, 
                       a_s = 2.5, b_s = 0.025,
                       seed = seed), 
          fixed_p = 5)
#z
z$acc/N
z$time
chain = unlist(z$chain)
sum( is.na(chain) )
chain[2, ] = tanh( chain[2, ] )
chain[3, ] = exp( chain[3, ] )
#mu = -1
#phi = 0.95
#sigma = 0.25
############################### Convergence analysis
# Trace plots
burn = 0
burned = chain[, -c(1:burn)]
par(mfrow = c(2, 2))
plot(burned[1, ], type = 'l')
plot(acf(burned[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[2, ], type = 'l')
plot(acf(burned[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[3, ], type = 'l')
plot(acf(burned[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

# Jumps
lags = 0
jumps = seq(1, N - burn, by = lags)
burned_lag = burned[, jumps]
par(mfrow = c(2, 2))
plot(burned_lag[1, ], type = 'l')
plot(acf(burned_lag[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[2, ], type = 'l')
plot(acf(burned_lag[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[3, ], type = 'l')
plot(acf(burned_lag[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

plot(burned[1, ], burned[2, ], xlab = 'mu', ylab = 'phi')
plot(burned[1, ], burned[3, ], xlab = 'mu', ylab = 'sigma')
plot(burned[2, ], burned[3, ], xlab = 'phi', ylab = 'sigma')
###########################################################################
N_new = length( burned )
mcmcchain = coda::as.mcmc( t(burned) )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff = coda::effectiveSize(mcmcchain)
IF = N_new / N_eff
IF
# |G| > 1.96 evidencia não convergencia
coda::geweke.diag(mcmcchain)

plot(burned_lag[1, ], burned_lag[2, ], xlab = 'mu', ylab = 'phi')
plot(burned_lag[1, ], burned_lag[3, ], xlab = 'mu', ylab = 'sigma')
plot(burned_lag[2, ], burned_lag[3, ], xlab = 'phi', ylab = 'sigma')

par(mfrow = c(1, 3))
hist(burned_lag[1, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[2, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[3, ], breaks = 30, xlab = '', main = '')
par(mfrow = c(1, 1))

data = matrix(c( mean(burned_lag[1, ]), 
                 quantile(burned_lag[1, ], probs = c(0.025, 0.975)),
                 mean(burned_lag[2, ]),
                 quantile(burned_lag[2, ], probs = c(0.025, 0.975)),
                 mean(burned_lag[3, ]),
                 quantile(burned_lag[3, ], probs = c(0.025, 0.975)) ), 
              nrow = 3, byrow = TRUE)
row.names(data) = c('mu_hat', 'phi_hat', 'sigma_hat')
colnames(data) = c('mean', '5%', '97,5%')
data










