################################################################################
#### librarys
library(ggplot2)
path = 'C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/normal/normal_model.cpp'
Rcpp::sourceCpp( path )

N = 3e4
samples = svmn(N,
               L_theta = 20, eps_theta = 0.5,
               L_b = 20, eps_b = 0.1,
               L_h = 50, eps_h = 0.015,
               y_T = c(0, y),
               seed = 976335 ) 

samples$time / 60
samples$acc / N
################## Save outputs
#save(samples, file = 'C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/normal/ES.RDara')
load('C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/normal/ES.RDara')
chain_theta = samples$chain$chain_theta
chain_b     = samples$chain$chain_b
chain_h     = samples$chain$chain_h
# Transformations
chain_theta[2, ] = tanh( chain_theta[2, ] )
chain_theta[3, ] = exp( chain_theta[3, ] )
chain_b[2, ]     = tanh( chain_b[2, ] )
############################### Convergence analysis
################### Trace plots
### burn
burn = 1e4
chain_theta  = chain_theta[, - c( 1:burn ) ] 
chain_b     = chain_b[, - c( 1:burn ) ]
chain_h     = chain_h[, - c( 1:burn ) ]
# Jumps
lags = 10
jumps = seq(1, N - burn, by = lags)
chain_theta = chain_theta[, jumps ]
chain_b     = chain_b[, jumps ]
chain_h     = chain_h[, jumps ]
N_new = length( jumps )
###############################################################################
###############################################################################
############################## Análise gráfica
mat = matrix(seq(1, 18), nrow = 3, ncol = 6)
layout( mat )
############################### theta
plot(chain_theta[1, ], type = 'l', main = 'mu', xlab = '', ylab = '')
plot(acf(chain_theta[1, ], lag.max = 100, plot = FALSE)[1:100], main = '', xlab = '', ylab = '')
hist(chain_theta[1, ], breaks = 40,  main = '', xlab = '', ylab = '')
plot(chain_theta[2, ], type = 'l', main = 'phi', xlab = '', ylab = '')
plot(acf(chain_theta[2, ], lag.max = 100, plot = FALSE)[1:100], main = '',xlab = '', ylab = '')
hist(chain_theta[2, ], breaks = 40, main = '',xlab = '', ylab = '')
plot(chain_theta[3, ], type = 'l', main = 'sigma', xlab = '', ylab = '')
plot(acf(chain_theta[3, ], lag.max = 100, plot = FALSE)[1:100], main = '',xlab = '', ylab = '')
hist(chain_theta[3, ], breaks = 40, main = '',xlab = '', ylab = '')
############################### b
plot(chain_b[1, ], type = 'l', main = 'b0', xlab = '', ylab = '')
plot(acf(chain_b[1, ], lag.max = 100, plot = FALSE)[1:100], main = '', xlab = '', ylab = '')
hist(chain_b[1, ], breaks = 40, main = '', xlab = '', ylab = '')
plot(chain_b[2, ], type = 'l', main = 'b1', xlab = '', ylab = '')
plot(acf(chain_b[2, ], lag.max = 100, plot = FALSE)[1:100], main = '', xlab = '', ylab = '')
hist(chain_b[2, ], breaks = 40, main = '', xlab = '', ylab = '')
plot(chain_b[3, ], type = 'l', main = 'b2', xlab = '', ylab = '')
plot(acf(chain_b[3, ], lag.max = 100, plot = FALSE)[1:100], main = '', xlab = '', ylab = '')
hist(chain_b[3, ], breaks = 40, main = '', xlab = '', ylab = '')
############################## Análise numérica
############################### theta
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
mc_error_theta = round( apply( chain_theta, 
                               MARGIN = 1, 
                               FUN = sd) / sqrt( N_eff_theta ), 
                        5 )
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
############################### b
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
mc_error_b = round( apply( chain_b, 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_b ), 
                    5 )
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
# Summary Table ####
data = data_theta
data = rbind( data, data_b )
#colnames( data ) = c('vdd', 'média', 'sd', '2.5%', '97.5%', 'CD', 'IF', 'mc_error')
colnames( data ) = c('média', 'sd', '2.5%', '97.5%', 'CD', 'IF', 'mc_error')
data = round( data, 4 )
data
###############################################################################
###############################################################################
############################### h
h_hat = apply(chain_h, MARGIN = 1, FUN = mean)
h_min = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.025) )
h_max = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:T, h, h_hat, h_min, h_max), ncol = 5)
#data = matrix(c(1:T, h_hat, h_min, h_max), ncol = 4)
data = data.frame(data)
names(data) = c('obs', 'vdd', 'media', 'min','max')
#names(data) = c('obs', 'media', 'min','max')
#plots
#a = sample(1:(T - 151), 1)
#g = ggplot(data[ a:(250 + a), ])
g = ggplot(data)
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
theta_hat = c( b_hat, h_hat )
theta_draws = chain_b
theta_draws = rbind( theta_draws, chain_h )
theta_draws = as.matrix( theta_draws )
############### dic deviance information criterion:
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# theta = b, h
# data = c( y0, y )
# theta_hat = ( b_hat, h_hat )
# theta_draws = burned_lag
log_lik = function(theta_t, data){
  # função checada (13/04/23)
  T = length( data ) - 1
  
  b0_t = theta_t[1]
  b1_t = theta_t[2]
  b2_t = theta_t[3]
  h_t = theta_t[4:(T+3)]
  
  log_l = dnorm(data[2:( T + 1 )], 
                mean = b0 + b1 * data[1:T] + b2 * exp( h_t ), 
                sd = exp( 0.5 * h_t ), 
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
dic_normal = dic(data = c(0, y), 
                 theta_draws = theta_draws, 
                 theta_hat )

############### loo
lik = function(data_i, draws, data_, data_past){
  #data_ = ( y_{1}, y_{2}, ..., y_{T} )
  #data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  k = which( data_ == as.numeric( data_i ) )
  log_l = NULL
  N = ncol( draws )
  
  for(col in 1:N){
    
    b0_draws = draws[1, col]
    b1_draws = draws[2, col]
    b2_draws = draws[3, col]
    h_draws  = draws[3 + k, col]
    
    log_l[col] = dnorm(data_i, mean = b0_draws + b1_draws * data_past[k] + b2_draws * exp( h_draws ), 
                       sd = exp( 0.5 * h_draws ) )
                       
  }
  
  return( log_l )
}

r_eff = loo::relative_eff(lik,
                          chain_id = rep(1, ncol( theta_draws ) ),
                          data = matrix( y , ncol = 1), 
                          draws = theta_draws,
                          data_ = y,
                          data_past = c( y0, y[1:(T-1)] )
                          #cores = getOption('mc.cores', 3)
)

# or set r_eff = NA
loo_normal = loo::loo(lik, 
                      #r_eff = NA,
                      r_eff = r_eff, 
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
  
  for(col in 1:N){
    
    b0_draws = draws[1, col]
    b1_draws = draws[2, col]
    b2_draws = draws[3, col]
    h_draws  = draws[3 + k, col]
    
    log_l[col] = dnorm(data_i, mean = b0_draws + b1_draws * data_past[k] + b2_draws * exp( h_draws ), 
                       sd = exp( 0.5 * h_draws ), log = TRUE )
    
  }
  
  return( log_l )
}

waic_normal = loo::waic(log_lik, 
                        data = matrix( y, ncol = 1 ), 
                        draws = theta_draws,
                        data_ = y,
                        data_past = c( y0, y[1:(T-1)] )
                        )

normal_criterium = matrix( c(dic_normal, 
                            loo_normal$looic,
                            waic_normal$waic), nrow = 1)
row.names( normal_criterium ) = c('normal')
colnames( normal_criterium ) = c('dic', 'loo', 'waic')
normal_criterium
