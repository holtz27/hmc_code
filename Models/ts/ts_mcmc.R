################################################################################
#### librarys
library( ggplot2 )
Rcpp::sourceCpp( 'ts_model.cpp' )

# Sampling
N = 5e4
samples = svm_smn_ts(N,
                     L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                     L_b = 20, eps_b = 0.1, 
                     L_h = 50, eps_h = 0.015,
                     L_v = 20, eps_v = 0.025, 
                     y_T = c(y0, y), 
                     seed = 0 )
samples$acc / N
samples$time / 60

################## Save outputs
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
burn = 2e3
# Jumps
lags = 10
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
par( mfrow = c(2, 3) )
plot(chain_theta[1, ], type = 'l', main = '', xlab = '', ylab = 'mu')
plot(acf(chain_theta[1, ], lag.max = 200, plot = FALSE)[1:200], main = '', 
     xlab = '', ylab = '')
hist(chain_theta[1, ],  main = '', xlab = '', ylab = '', breaks = 40)
plot(chain_theta[2, ], type = 'l', main = '', xlab = '', ylab = 'phi')
plot(acf(chain_theta[2, ], lag.max = 200, plot = FALSE)[1:200], main = '', 
     xlab = '', ylab = '')
hist(chain_theta[2, ], main = '', xlab = '', ylab = '', breaks = 40)
plot(chain_theta[3, ], type = 'l',  main = '', xlab = '', ylab = 'sigma')
plot(acf(chain_theta[3, ], lag.max = 200, plot = FALSE)[1:200], main = '', 
     xlab = '', ylab = '')
hist(chain_theta[3, ], main = '', xlab = '', ylab = '', breaks = 40)
par( mfrow = c(1, 1) )

############### Numeric Analysis
mcmcchain_theta = coda::as.mcmc( t( chain_theta ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmcchain_theta )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmcchain_theta )
IF_theta = N_new / N_eff_theta
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_theta = apply( chain_theta, 
                        MARGIN = 1, 
                        FUN = sd) / sqrt( N_eff_theta )
                       
theta_hat = apply( chain_theta, MARGIN = 1, FUN = mean )
theta_sd = apply( chain_theta, MARGIN = 1, FUN = sd )
theta_min = apply( chain_theta, MARGIN = 1, FUN = quantile, probs = c(0.025) )
theta_max = apply( chain_theta, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data_theta = matrix(
                    c(theta_hat,
                      theta_sd,
                      theta_min,
                      theta_max,
                      CD_theta$z,
                      IF_theta,
                      mc_error_theta), nrow = 3, byrow = FALSE
                    )
row.names( data_theta ) = c('mu', 'phi', 'sigma')
###############################################################################
###############################################################################
############################### b
par( mfrow = c(2, 3) )
plot(chain_b[1, ], type = 'l', main = '', xlab = '', ylab = 'b0')
plot(acf(chain_b[1, ], lag.max = 100, plot = FALSE)[1:100], main = '', 
     xlab = '', ylab = '')
hist(chain_b[1, ], main = '', xlab = '', ylab = '', breaks = 40)
plot(chain_b[2, ], type = 'l', main = '', xlab = '', ylab = 'b1')
plot(acf(chain_b[2, ], lag.max = 100, plot = FALSE)[1:100], main = '', 
     xlab = '', ylab = '')
hist(chain_b[2, ], main = '', xlab = '', ylab = '', breaks = 40)
plot(chain_b[3, ], type = 'l', main = '', xlab = '', ylab = 'b2')
plot(acf(chain_b[3, ], lag.max = 100, plot = FALSE)[1:100], main = '', 
     xlab = '', ylab = '')
hist(chain_b[3, ], main = '', xlab = '', ylab = '', breaks = 40)
par( mfrow = c(1, 1) )

############### Numeric Analysis
mcmcchain_b = coda::as.mcmc( t( chain_b ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_b = coda::geweke.diag( mcmcchain_b )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_b = coda::effectiveSize( mcmcchain_b )
IF_b = N_new / N_eff_b
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_b = apply( chain_b, 
                    MARGIN = 1, 
                    FUN = sd) / sqrt( N_eff_b )

b_hat = apply( chain_b, MARGIN = 1, FUN = mean )
b_sd = apply( chain_b, MARGIN = 1, FUN = sd )
b_min = apply( chain_b, MARGIN = 1, FUN = quantile, probs = c(0.025) )
b_max = apply( chain_b, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data_b = matrix(
                c(b_hat,
                  b_sd,
                  b_min,
                  b_max,
                  CD_b$z,
                  IF_b,
                  mc_error_b), nrow = 3, byrow = FALSE
               )
row.names( data_b ) = c('b0', 'b1', 'b2')
###############################################################################
###############################################################################
############################### e = log( v )
par( mfrow = c(1, 3) )
plot(chain_e, type = 'l', main = '', xlab = '', ylab = 'e = log( v )')
#abline( h = log( v ) )
plot(acf(chain_e, lag.max = 100, plot = FALSE)[1:100], main = '', 
     xlab = '', ylab = '')
hist(chain_e, main = '', xlab = '', ylab = '', breaks = 40)
par( mfrow = c(1, 1) )

############### Numeric Analysis
mcmcchain_e = coda::as.mcmc( chain_e ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_e = coda::geweke.diag( mcmcchain_e )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_e = coda::effectiveSize( mcmcchain_e )
IF_e = N_new / N_eff_e
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_e = sd( chain_e ) / sqrt( N_eff_e )
e_hat = mean( chain_e )
e_sd = sd( chain_e )
e_min = quantile( chain_e, probs = c(0.025) )
e_max = quantile( chain_e, probs = c(0.975) )
data_e = matrix(
                c(e_hat,
                  e_sd, 
                e_min,
                e_max,
                CD_e$z,
                IF_e,
                mc_error_e), nrow = 1, byrow = FALSE
               )
row.names( data_e ) = c('e')
###############################################################################
###############################################################################
# Summary Table
data = data_theta
data = rbind( data, data_b, data_e )
#data = cbind( c(mu, phi, sigma, b0, b1, b2, log(v) ), data )
#colnames( data ) = c('vdd', 'média', 'sd', '2.5%', '97.5%', 'CD', 'IF', 'mc_error')
colnames( data ) = c('média', 'sd', '2.5%', '97.5%', 'CD', 'IF', 'mc_error')
data = round( data, 4 )
data
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
a = sample(1:(T - 101), 1)
f = ggplot(data2[a:(a + 100), ]) 
f = f + geom_line(aes(obs, media))
f = f + geom_line(aes(obs, vdd), color = 'red')
f = f + geom_line(aes(obs, min), linetype = 'dashed')
f = f + geom_line(aes(obs, max), linetype = 'dashed')
f
############### Numeric Analysis
mcmcchain_l = coda::as.mcmc( t( chain_l ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_l = coda::geweke.diag( mcmcchain_l )
# Fração de valores que est]ao no intervalo ( -1.96 , 1.96 )
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
#plots
a = sample(1:(T - 251), 1)
g = ggplot(data[ a:(250 + a), ]) 
g = g + geom_line(aes(obs, media))
g = g + geom_line(aes(obs, vdd), color = 'red')
g = g + geom_line(aes(obs, min), linetype = 'dashed')
g = g + geom_line(aes(obs, max), linetype = 'dashed')
g
############### Numeric Analysis
mcmcchain_h = coda::as.mcmc( t( chain_h ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_h = coda::geweke.diag( mcmcchain_h )
# Fração de valores que est]ao no intervalo ( -1.96 , 1.96 )
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
###############################################################################
###############################################################################
############################## Model Selection
# construindo theta_hat e theta_draws
theta_hat = c( b_hat, h_hat, l_hat )
theta_draws = chain_b
theta_draws = rbind( theta_draws, chain_h )
theta_draws = rbind( theta_draws, chain_l )
############### dic deviance information criterion:
# p( y | theta ) = p( y | b, theta_h, v, h, l) = p( y | b, h, l )
# theta = b, h, l
# data = c( y0, y )
# theta_hat = ( b_hat, h_hat, l_hat )
# theta_draws = burned_lag
log_lik = function(theta_t, data){
  # função checada (13/04/23)
  T = length( data ) - 1
  
  b0_t = theta_t[1]
  b1_t = theta_t[2]
  b2_t = theta_t[3]
  h_t = theta_t[4:(T+3)]
  l_t = theta_t[(T + 4):(2 * T + 3)]
  log_l = dnorm(data[2:( T + 1 )], 
                mean = b0 + b1 * data[1:T] + b2 * exp( h_t ), 
                sd = exp( 0.5 * h_t ) / sqrt( l_t ), 
                log = TRUE)
  
  return( sum( log_l ) )
  
}
log_lik_i = function(data, theta_draws){
  # função checada (13/04/23)
  x = apply(X = theta_draws, MARGIN = 2, FUN = log_lik, data)
  
  return( x )
}
dic = function(data, theta_draws, theta_hat){
  
  pD = 2 * ( log_lik( theta_hat, data ) - mean( log_lik_i(data, theta_draws) ) )      
  DIC = - 2 * log_lik( theta_hat, data ) + 2 * pD  
  
  return( DIC )
} 
# calculando DIC
svm_ts = dic( data = c(y0, y), 
              theta_draws = theta_draws, 
              theta_hat )
############### loo
lik = function(data_i, draws, data_, data_past){
  #data_ = ( y_{1}, y_{2}, ..., y_{T} )
  #data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  k = which( data_ == as.numeric( data_i ) )
  log_l = NULL
  N = ncol( draws )
  T = 0.5 * (nrow( draws ) - 3 )
  for(col in 1:N){
    
    b0_draws = draws[1, col]
    b1_draws = draws[2, col]
    b2_draws = draws[3, col]
    h_draws  = draws[3 + k, col]
    l_draws  = draws[T + 3 + k, col] 
    
    log_l[col] = dnorm(data_i, mean = b0_draws + b1_draws * data_past[k] + b2_draws * exp( h_draws ), 
                       sd = exp( 0.5 * h_draws ) / sqrt( l_draws ) )
    
  }
  
  return( log_l )
}

r_eff = loo::relative_eff(lik,
                          chain_id = rep(1, ncol( theta_draws ) ),
                          data = as.matrix( y ), 
                          draws = theta_draws,
                          data_ = y,
                          data_past = c( y0, y[1:(T-1)] ),
                          cores = getOption('mc.cores', 3)
)

# or set r_eff = NA
loo::loo(lik, 
         r_eff = NA,
         #r_eff = r_eff, 
         data = as.matrix( y ), 
         draws = theta_draws,
         data_ = y,
         data_past = c( y0, y[1:(T-1)] ),
         cores = getOption('mc.cores', 3)
)
############### waic
log_lik = function(data_i, draws, data_, data_past){
  #data_ = ( y_{1}, y_{2}, ..., y_{T} )
  #data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  k = which( data_ == as.numeric( data_i ) )
  log_l = NULL
  N = ncol( draws )
  T = 0.5 * (nrow( draws ) - 3 )
  for(col in 1:N){
    
    b0_draws = draws[1, col]
    b1_draws = draws[2, col]
    b2_draws = draws[3, col]
    h_draws  = draws[3 + k, col]
    l_draws  = draws[T + 3 + k, col] 
    
    log_l[col] = dnorm(data_i, mean = b0_draws + b1_draws * data_past[k] + b2_draws * exp( h_draws ), 
                       sd = exp( 0.5 * h_draws ) / sqrt( l_draws ), log = TRUE )
    
  }
  
  return( log_l )
}

loo::waic(log_lik, 
          data = matrix( y, ncol = 1 ), 
          draws = theta_draws,
          data_ = y,
          data_past = c( y0, y[1:(T-1)] )
)
