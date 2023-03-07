path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/Parâmetros/b/jacobiano/b_jacobiano.cpp'
Rcpp::sourceCpp(path)

N = 1
z = rmhmc(eps = c(0.05, 0.1, 0.05), 
          min_L = 40, max_L = 50, 
          theta_init = c(mean(y), atanh(0.8), 2), 
          param = list(y_T = c(y0, y), 
                       h = rnorm(T), 
                       l = rgamma(T, shape = 2), 
                       mu_b0 = 0, s_b0 = 10, 
                       a_b1 = 5, b_b1 = 1.5, 
                       mu_b2 = 0, s_b2 = 10,
                       seed = 42), 
          fixed_p = 5)
z
z$acc/N
z$time

chain = unlist(z$chain)
chain[2, ] = tanh(chain[2, ])

#b0 = 0.001, b1 = 0.2, b2 = -0.2
############################### Convergence analysis
############### Análise gráfica
burn = 10
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
lags = 3
jumps = seq(1, N - burn, by = 1 + lags)
burned_lag = burned[, jumps]

par(mfrow = c(2, 2))
plot(burned_lag[1, ], type = 'l')
plot(acf(burned_lag[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[2, ], type = 'l')
plot(acf(burned_lag[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[3, ], type = 'l')
plot(acf(burned_lag[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

############### Análise numérica
mcmcchain = coda::as.mcmc( t( burned_lag ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD = coda::geweke.diag(mcmcchain)
CD
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_new = length( burned_lag[1, ] )
N_eff = coda::effectiveSize(mcmcchain)
IF = N_new / N_eff
IF
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error = round( apply( burned_lag, MARGIN = 1, FUN = sd) / sqrt( N_eff ), 5 )
mc_error

# Histogramas
par(mfrow = c(1, 3))
hist(burned_lag[1, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[2, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[3, ], breaks = 30, xlab = '', main = '')
par(mfrow = c(1, 1))
data = matrix(c( mean(burned_lag[1, ]), 
                 quantile(burned_lag[1, ], probs = c(0.05, 0.975)),
                 mean(burned_lag[2, ]),
                 quantile(burned_lag[2, ], probs = c(0.05, 0.975)),
                 mean(burned_lag[3, ]),
                 quantile(burned_lag[3, ], probs = c(0.05, 0.975)) ), 
              nrow = 3, byrow = TRUE)
row.names(data) = c('b0_hat', 'b1_hat', 'b2_hat')
colnames(data) = c('mean', '5%', '97,5%')
data










