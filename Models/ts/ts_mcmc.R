
library( ggplot2 )

Rcpp::sourceCpp('ts_model.cpp')

# Sampling
N = 5e3
samples = svm_smn_ts(N,
                     L_theta = 100, eps_theta = 0.5, 
                     L_b = 20, eps_b = 0.1, 
                     L_h = 50, eps_h = 0.015,
                     L_v = 80, eps_v = 0.025, 
                     y_T = c(y0, y), 
                     seed = 0 )

samples$acc / N
samples$time / 60

# Saving
#save(samples, file = 'ts_samples.RDara')
#load('ts_samples.RDara')

chain_theta = samples$chain$chain_theta
chain_b = samples$chain$chain_b
chain_h = samples$chain$chain_h
chain_e = samples$chain$chain_v
chain_l = samples$chain$chain_l

# Transformations
chain_theta[2, ] = tanh( chain_theta[2, ] )
chain_theta[3, ] = exp( chain_theta[3, ] )
chain_b[2, ]     = tanh( chain_b[2, ] )
chain_v          = exp( chain_e[1, ] )

############################### Convergence analysis
################### Trace plots
### burn
burn = 0
# Jumps
lags = 1
jumps = seq(1, N - burn, by = lags)

chain_theta  = chain_theta[, - c( 1:burn ) ] 
chain_b     = chain_b[, - c( 1:burn ) ]
chain_h     = chain_h[, - c( 1:burn ) ]
chain_e     = chain_e[ - c( 1:burn ) ]
chain_v     = chain_v[ - c( 1:burn ) ]
chain_l     = chain_l[, - c( 1:burn ) ]  

chain_theta = chain_theta[, jumps ]
chain_b     = chain_b[, jumps ]
chain_h     = chain_h[, jumps ]
chain_e     = chain_e[ jumps ]
chain_v     = chain_v[ jumps ]
chain_l     = chain_l[, jumps ]

N_new = length( jumps )
###############################################################################
###############################################################################
############################### theta

par(mfrow = c(2, 3))
plot(chain_theta[1, ], type = 'l', main = 'mu', ylab = '')
plot(acf(chain_theta[1, ], lag.max = 100, plot = FALSE)[1:100])
hist(chain_theta[1, ], breaks = 40)
plot(chain_theta[2, ], type = 'l', main = 'phi', ylab = '')
plot(acf(chain_theta[2, ], lag.max = 100, plot = FALSE)[1:100])
hist(chain_theta[2, ], breaks = 40)
plot(chain_theta[3, ], type = 'l', main = 'sigma', ylab = '')
plot(acf(chain_theta[3, ], lag.max = 100, plot = FALSE)[1:100])
hist(chain_theta[3, ], breaks = 40)
par(mfrow = c(1, 1))

############### Análise numérica
mcmcchain_theta = coda::as.mcmc( t( chain_theta ) )
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
mc_error_theta = round( apply( chain_theta, 
                               MARGIN = 1, 
                               FUN = sd) / sqrt( N_eff_theta ), 
                        5 )
mc_error_theta

###############################################################################
###############################################################################
############################### b

par( mfrow = c(2, 3) )
plot(chain_b[1, ], type = 'l', main = 'b0', ylab = '')
plot(acf(chain_b[1, ], lag.max = 100, plot = FALSE)[1:100])
hist(chain_b[1, ], breaks = 40)
plot(chain_b[2, ], type = 'l', main = 'b1', ylab = '')
plot(acf(chain_b[2, ], lag.max = 100, plot = FALSE)[1:100])
hist(chain_b[1, ], breaks = 40)
plot(chain_b[3, ], type = 'l', main = 'b2', ylab = '')
plot(acf(chain_b[3, ], lag.max = 100, plot = FALSE)[1:100])
hist(chain_b[1, ], breaks = 40)
par( mfrow = c(1, 1) )

############### Análise numérica
mcmcchain_b = coda::as.mcmc( t( chain_b ) )
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
mc_error_b = round( apply( chain_b, 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_b ), 
                    5 )
mc_error_b

###############################################################################
###############################################################################
############################### v

par(mfrow = c(1, 3))
plot(chain_e, type = 'l', main = 'e = log( v )')
abline( h = log( v ) )
plot(acf(chain_e, lag.max = 100, plot = FALSE)[1:100])
hist(chain_e, breaks = 40)
par(mfrow = c(1, 1))

############### Análise numérica
mcmcchain_e = coda::as.mcmc( chain_e ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_e = coda::geweke.diag( mcmcchain_e )
CD_e
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_e = coda::effectiveSize( mcmcchain_e )
IF_e = N_new / N_eff_e
IF_e
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_e = round( sd( chain_e ) / sqrt( N_eff_e ), 5 )
mc_error_e

###############################################################################
###############################################################################
############################### l
 
l_hat = apply(chain_l, MARGIN = 1, FUN = mean)
l_min = apply(chain_l, MARGIN = 1, FUN = quantile, probs = c(0.025) )
l_max = apply(chain_l, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data2 = matrix(c(1:T, l, l_hat, l_min, l_max), ncol = 5)
data2 = data.frame(data2)
names(data2) = c('obs', 'vdd', 'media', 'min', 'max')
# plot1

a = sample(1:T, 1)
f = ggplot(data2[a:(a + 100), ]) 
f = f + geom_line(aes(obs, media))
f = f + geom_line(aes(obs, vdd), color = 'red')
f = f + geom_line(aes(obs, min), linetype = 'dashed')
f = f + geom_line(aes(obs, max), linetype = 'dashed')
f

############### Análise numérica
mcmcchain_l = coda::as.mcmc( t( chain_l ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_l = coda::geweke.diag( mcmcchain_l )
# Espera-se que mais de 95% desdes valores estejam no intervalo ( -1.96 , 1.96 )
geweke_l = sum( abs( CD_l$z ) < 1.96 ) / T
geweke_l
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_l = coda::effectiveSize( mcmcchain_l )
IF_l = N_new / N_eff_l
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_l = round( apply( chain_l, 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_l ), 5 )
# plots
par( mfrow = c(1,3) )
plot( CD_l$z, main = 'Geweke diagnostic', xlab = '', ylab = '' )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_l, main = 'Inefficiency factors', xlab = '', ylab = '' )
abline(h = 1)
plot( mc_error_l, main = 'MCMC errors', xlab = '', ylab = '' )
par( mfrow = c(1,1) )

###############################################################################
###############################################################################
############################### h

h_hat = apply(chain_h, MARGIN = 1, FUN = mean)
h_min = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.025) )
h_max = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:T, h, h_hat, h_min, h_max), ncol = 5)
data = data.frame(data)
names(data) = c('obs', 'vdd', 'media', 'min','max')

a = sample(1:T, 1)

g = ggplot(data[ a:(250 + a), ]) 
g = g + geom_line(aes(obs, media))
g = g + geom_line(aes(obs, vdd), color = 'red')
g = g + geom_line(aes(obs, min), linetype = 'dashed')
g = g + geom_line(aes(obs, max), linetype = 'dashed')
g

############### Análise numérica
mcmcchain_h = coda::as.mcmc( t( chain_h ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_h = coda::geweke.diag( mcmcchain_h )
# Espera-se que 95% desdes valores estejam no intervalo ( -1.96 , 1.96 )
geweke_h = sum( abs( CD_h$z ) < 1.96 ) / T
geweke_h
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_h = coda::effectiveSize( mcmcchain_h )
IF_h = N_new / N_eff_h
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_h = round( apply( chain_h, 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_h ), 
                    5 )

# plots
par( mfrow = c(1,3) )
plot( CD_h$z, main = 'Geweke diagnostic', xlab = '', ylab = '' )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_h, main = 'Inefficiency factors', xlab = '', ylab = '' )
abline(h = 1)
plot( mc_error_h, main = 'MCMC errors', xlab = '', ylab = '' )
par( mfrow = c(1,1) )
