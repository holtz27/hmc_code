################################################################################
#### librarys
library(ggplot2)
################################################################################
path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/Modelos/ts/jacobiano/ts_partial_model.cpp'
Rcpp::sourceCpp(path)
################################################################################
y_T = c(y0, y)
theta_init = c( 1, atanh( 0.86 ), log( 0.5 ) )
b_init = c( mean(y), atanh( 0.8 ), 2 )
e_init = log( 40 )

init = c(theta_init,
         b_init,
         e_init)

param = list(y_T = c(y0, y),
             theta = theta_init,
             b = b_init,
             h = h, 
             l = l,
             e = e_init, 
             mu_0 = 0, s_0 = 3.2, 
             a_phi = 20, b_phi = 1.5, 
             a_s = 2, b_s = 10,
             mu_b0 = 0, s_b0 = 3.2, 
             a_b1 = 5, b_b1 = 1.5, 
             mu_b2 = 0, s_b2 = 3.2,
             a_v = 1, b_v = 0.1,
             seed = 0)
N = 1e3
x = svm_smn(N, 
            eps_theta = c(0.1, 0.1, 0.1), min_L_theta = 20, max_L_theta = 20,
            eps_b = c(0.1, 0.1, 0.1), min_L_b = 20, max_L_b = 20,
            eps_e = 0.1, min_L_e = 50, max_L_e = 50,
            init = init, 
            param)

x$acc / N
#x$time
chain = unlist(x$chain)
sum( is.na( chain ) )

theta_init = chain[1:3, N + 1]
b_init = chain[4:6, N + 1]
e_init = chain[7, N + 1]
init = c(theta_init,
         b_init,
         e_init)

# Transformations
chain[2, ] = tanh( chain[2, ] )
chain[3, ] = exp( chain[3, ] )
chain[5, ] = tanh( chain[5, ] )
chain[7, ] = exp( chain[7, ] )
############################### Convergence analysis
################### Trace plots
### burn
burn = 0
burned = as.matrix(chain[, -c(1:burn)])
###############################################################################
###############################################################################
############################### theta

par(mfrow = c(2, 2))
plot(burned[1, ], type = 'l', main = 'mu')
plot(acf(burned[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[2, ], type = 'l', main = 'phi')
plot(acf(burned[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[3, ], type = 'l', main = 'sigma')
plot(acf(burned[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))
# Jumps
lags = 1
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
mc_error_theta = round( apply( burned_lag[1:3, ], MARGIN = 1, FUN = sd) / sqrt( N_eff_theta ), 5 )
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
mc_error_b = round( apply( burned_lag[4:6, ], MARGIN = 1, FUN = sd) / sqrt( N_eff_b ), 5 )
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
############################### h
X = burned_lag[7:(T + 6), ]
h_hat = apply(X, MARGIN = 1, FUN = mean)
h_min = apply(X, MARGIN = 1, FUN = quantile, probs = c(0.025) )
h_max = apply(X, MARGIN = 1, FUN = quantile, probs = c(0.975) )
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
mcmcchain_h = coda::as.mcmc( t( X ) ) 
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
mc_error_h = round( apply( burned_lag[7:(T + 6), ], MARGIN = 1, FUN = sd) / sqrt( N_eff_h ), 5 )
plot( mc_error_h )

###############################################################################
###############################################################################
############################### l
Y = burned_lag[(T + 7):(2*T + 6), ]
l_hat = apply(Y, MARGIN = 1, FUN = mean)
l_min = apply(Y, MARGIN = 1, FUN = quantile, probs = c(0.025) )
l_max = apply(Y, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data2 = matrix(c(1:T, l, l_hat, l_min, l_max), ncol = 5)
data2 = data.frame(data2)
names(data2) = c('obs', 'vdd', 'media', 'min', 'max')
# plot1
f = ggplot(data2) 
f = f + geom_line(aes(obs, media))
f = f + geom_line(aes(obs, vdd), color = 'red')
f = f + geom_line(aes(obs, min), linetype = 'dashed')
f = f + geom_line(aes(obs, max), linetype = 'dashed')
f

############### Análise numérica
mcmcchain_l = coda::as.mcmc( t( Y ) ) 
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
mc_error_l = round( apply( Y, MARGIN = 1, FUN = sd) / sqrt( N_eff_l ), 5 )
plot( mc_error_l )




