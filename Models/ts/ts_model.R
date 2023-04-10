

Rcpp::sourceCpp('run_svm_smn_ts.cpp')

theta_init = c(0.005,
               0.5 * ( log( 1 + 0.98 ) - log( 1 - 0.98 ) ),
               log( sqrt( 0.017 ) ) )

b_init = c(0.3,
           0.5 * ( log( 1 + 0.03 ) - log( 1 - 0.03 ) ),
           -0.025 )

h_init = matrix( nrow = T, ncol = 1 )
h_init[1, 1] = 0.005 + sqrt( 0.03 ) / (1 - 0.95 * 0.95 ) * rnorm(1)
for( kt in 2:T ){
  h_init[ kt, 1 ] = 0.005 + 0.95 * ( h_init[ kt - 1, 1 ] - 0.005 ) + sqrt( 0.03 ) * rnorm(1);
}

v_init = log( 10 )

N = 5e4
samples = svm_smn_ts(N,
                     L_theta = 10, eps_theta = 0.01, theta_init,
                     L_b = 10, eps_b = 0.01, b_init,
                     L_h = 10, eps_h = 0.01, h_init,
                     L_v = 10, eps_v = 0.01, v_init, 
                     y_T = c(y0, y), 
                     seed = 1029384756 )

samples$acc / N
samples$time / 60

theta_init = samples$chain$chain_theta[, N + 1 ]
b_init     = samples$chain$chain_b[, N + 1 ]
h_init     = samples$chain$chain_h[, N + 1 ]
v_init     = samples$chain$chain_v[ N + 1 ]

chain_theta = samples$chain$chain_theta
chain_b = samples$chain$chain_b
chain_h = samples$chain$chain_h
chain_v = samples$chain$chain_v

# Transformations
chain_theta[2, ] = tanh( chain_theta[2, ] )
chain_theta[3, ] = exp( chain_theta[3, ] )
chain_b[2, ]     = tanh( chain_b[2, ] )
chain_v          = exp( chain_v )
############################### Convergence analysis
################### Trace plots
### burn
burn = 1e4

theta_burned = chain_theta[, - c( 1:burn ) ] 
chain_b      = chain_b[, - c( 1:burn ) ]
chain_h      = chain_h[, - c( 1:burn ) ]
chain_v      = chain_v[, - c( 1:burn ) ]
###############################################################################
###############################################################################
############################### theta

par(mfrow = c(2, 2))
plot(theta_burned[1, ], type = 'l', main = 'mu')
plot(acf(theta_burned[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(theta_burned[2, ], type = 'l', main = 'phi')
plot(acf(theta_burned[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(theta_burned[3, ], type = 'l', main = 'sigma')
plot(acf(theta_burned[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))
# Jumps
lags = 10
jumps = seq(1, N - burn, by = lags)
theta_burned_lag = burned[, jumps]
par(mfrow = c(2, 2))
plot(theta_burned_lag[1, ], type = 'l')
plot(acf(theta_burned_lag[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(theta_burned_lag[2, ], type = 'l')
plot(acf(theta_burned_lag[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(theta_burned_lag[3, ], type = 'l')
plot(acf(theta_burned_lag[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

N_new = length( burned_lag[1, ] )

############### Análise numérica
mcmcchain_theta = coda::as.mcmc( t( burned_lag[1:3, ] ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmcchain_theta )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmcchain_theta )
IF_theta = N_new / N_eff_theta
IF_theta
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_theta = round( apply( burned_lag[1:3, ], 
                               MARGIN = 1, 
                               FUN = sd) / sqrt( N_eff_theta ), 5 )
mc_error_theta

###############################################################################
###############################################################################
############################### b

par(mfrow = c(2, 2))
plot(burned_lag[4, ], type = 'l', main = 'b0')
plot(acf(burned_lag[4, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[5, ], type = 'l', main = 'b1')
plot(acf(burned_lag[5, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[6, ], type = 'l', main = 'b2')
plot(acf(burned_lag[6, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

############### Análise numérica
mcmcchain_b = coda::as.mcmc( t( burned_lag[4:6, ] ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_b = coda::geweke.diag( mcmcchain_b )
CD_b
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_b = coda::effectiveSize( mcmcchain_b )
IF_b = N_new / N_eff_b
IF_b
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_b = round( apply( burned_lag[4:6, ], 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_b ), 5 )
mc_error_b

###############################################################################
###############################################################################
############################### v

par(mfrow = c(1, 2))
plot(burned_lag[7, ], type = 'l', main = 'v')
plot(acf(burned_lag[7, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

############### Análise numérica
mcmcchain_v = coda::as.mcmc(  burned_lag[7, ] ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_v = coda::geweke.diag( mcmcchain_v )
CD_v
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_v = coda::effectiveSize( mcmcchain_v )
IF_v = N_new / N_eff_v
IF_v
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_v = round( sd( burned_lag[( 7), ] ) / sqrt( N_eff_v ), 5 )
mc_error_v

###############################################################################
###############################################################################
############################### l
L = burned_lag[8:(T + 7), ]
l_hat = apply(L, MARGIN = 1, FUN = mean)
l_min = apply(L, MARGIN = 1, FUN = quantile, probs = c(0.025) )
l_max = apply(L, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data2 = matrix(c(1:T, l, l_hat, l_min, l_max), ncol = 5)
data2 = data.frame(data2)
names(data2) = c('obs', 'vdd', 'media', 'min', 'max')
# plot1
f = ggplot(data2[151:250, ]) 
f = f + geom_line(aes(obs, media))
f = f + geom_line(aes(obs, vdd), color = 'red')
f = f + geom_line(aes(obs, min), linetype = 'dashed')
f = f + geom_line(aes(obs, max), linetype = 'dashed')
f

############### Análise numérica
mcmcchain_l = coda::as.mcmc( t( L ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_l = coda::geweke.diag( mcmcchain_l )
plot( CD_l$z )
abline(h = -1.96)
abline(h = 1.96)
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_l = coda::effectiveSize( mcmcchain_l )
IF_l = N_new / N_eff_l
plot( IF_l )
abline(h = 1)
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_l = round( apply( L, MARGIN = 1, FUN = sd) / sqrt( N_eff_l ), 5 )
plot( mc_error_l )


###############################################################################
###############################################################################
############################### h
H = burned_lag[(T + 8): (2 * T + 7), ]
h_hat = apply(H, MARGIN = 1, FUN = mean)
h_min = apply(H, MARGIN = 1, FUN = quantile, probs = c(0.025) )
h_max = apply(H, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:T, h, h_hat, h_min, h_max), ncol = 5)
data = data.frame(data)
names(data) = c('obs', 'vdd', 'media', 'min','max')
g = ggplot(data) 
g = g + geom_line(aes(obs, media))
g = g + geom_line(aes(obs, vdd), color = 'red')
g = g + geom_line(aes(obs, min), linetype = 'dashed')
g = g + geom_line(aes(obs, max), linetype = 'dashed')
g

############### Análise numérica
mcmcchain_h = coda::as.mcmc( t( H ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_h = coda::geweke.diag( mcmcchain_h )
plot( CD_h$z )
abline(h = -1.96)
abline(h = 1.96)
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_h = coda::effectiveSize( mcmcchain_h )
IF_h = N_new / N_eff_h
plot( IF_h )
abline(h = 1)
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_h = round( apply( burned_lag[(T + 8): (2 * T + 7), ], MARGIN = 1, FUN = sd) / sqrt( N_eff_h ), 5 )
plot( mc_error_h )
