###################################################################
########################## Modelo
# prioris: mu ~ N(m_0, s_0), sigma2 ~ Gama(a_0, b_0) 
# theta = (mu, sigma2)
############## Transformation sigma = exp(w) ######################
# theta = (mu, w)
###################################################################
###################################################################
logpost_theta = function(theta, param){
  # theta = (mu, w)
  
  mu = theta[1]
  w = theta[2]
  
  L = sum( dnorm(param$y, mu, exp(w), log = TRUE) )
  L = L + dnorm(mu, param$mu_0, param$s_0, log = TRUE)
  L = L + dgamma( exp( 2 * w ), shape = param$a, rate = param$b, log = TRUE)
  #x = x + 2*w
  
  return( L )
}
glogpost_theta = function(theta, param){
  mu = theta[1]
  w = theta[2]
  T = length(param$y)
  
  grad_mu = sum(param$y - mu) * exp( -2 * w )
  grad_mu = grad_mu - (mu - param$mu_0) / ( param$s_0 )**2
  
  grad_w = - T + sum( (param$y - mu)**2 ) * exp( -2 * w )
  grad_w = grad_w + 2 * (param$a - 1) - 2 * param$b * exp( 2 * w )
  
  return( matrix(c(grad_mu, grad_w), nrow = 2) )
}
G_theta_t = function(theta, param){
  mu = theta[1]
  w = theta[2]
  T = length(param$y)
  
  g11 =  T * exp(- 2 * w ) + 1 / param$s_0**2
  g22 = 2 * T  + 4 * param$b * exp(2 * w) 
  
  return( diag( c(g11, g22), 2, 2 ) )
}
d_G_mu = function(theta, param) return( matrix(0, 2, 2) )
d_G_sigma2 = function(theta, param){
  mu = theta[1]
  w = theta[2]
  T = length(param$y)
  
  g11 = -2 * T * exp( -2 * w )
  g22 = 8 * param$b * exp( 2 * w )
  
  return( diag(c(g11, g22), 2, 2) )
}
log_J_theta = function(theta, param){
  w = theta[2]
  return( log(2) + 2 * w )
}
###################################################################
########################## Algoritmo
###################################################################
H = function(p, theta, param){
  
  logpost = param$logpost
  G = param$G
  log_J = param$log_J
  
  p = matrix(p, ncol = 1)
  
  E = - logpost(theta, param) - log_J(theta, param)
  E = E + 0.5 * log( det( G(theta, param) ) )  
  E = E + 0.5 * t(p) %*% solve( G(theta, param) ) %*% p 
  
  return( E )
}
d_H_theta = function(p, theta, param){
  # param = (d_G = (d_G_theta1, ..., d_G_thetaN), ... )
  p = matrix(p, ncol = 1)
  D = length(p)
  d = matrix(nrow = D, ncol = 1)
  glogpost = param$glogpost
  G = param$G
  G.i = solve(G(theta, param))
  
  u = matrix(ncol = 1)
  
  for(i in 1:D){
    M = G.i %*% param$d_G[[i]](theta, param)
    u[i] = 0.5 * ( sum(diag( M )) - t(p) %*% M %*% d_H_p(p, theta, param) )  
    
  }
  
  u = matrix(u, ncol = 1)
  
  #d_H_theta
  d = - glogpost(theta, param) + u 
  
  return( d )
}
d_H_p = function(p, theta, param){
  p = matrix(p, ncol = 1)
  G = param$G
  return( solve( G(theta, param) ) %*% p )
}
glf = function(eps, L, theta_current, p_current, fixed_p, param){
  
  p_n = matrix(p_current, ncol = 1)
  theta_n = theta_hat = matrix(theta_current, ncol = 1)
  Proposal = matrix(nrow = nrow(p_n), ncol = 2)
  eps = matrix(eps, ncol = 1)
  
  for(t in 1:L){
    
    p_hat = p_n
    
    for(i in 1:fixed_p){
      p_hat = p_n - 0.5 * eps * d_H_theta(p_hat, theta_n, param)
    } 
    
    for(i in 1:fixed_p){
      theta_hat = theta_n + 0.5 * eps * ( d_H_p(p_hat, theta_n, param) + 
                                            d_H_p(p_hat, theta_hat, param) )
    }
    
    theta_n = theta_hat
    p_n = p_hat - 0.5 * eps * d_H_theta(p_hat, theta_n, param)
    
  }
  
  Proposal[, 1] = theta_n
  Proposal[, 2] = p_n
  
  return(Proposal)
}
rmhmc_R = function(N, eps, min_L, max_L, theta_init, param, fixed_p, 
                   seed = NULL){
  
  if ( !is.null(seed) ) set.seed(seed)
  time.init = Sys.time()
  acc = 0;
  theta_current = matrix(theta_init, ncol = 1)
  r = nrow(theta_current)
  logpost = param$logpost
  glogpost = param$glogpost
  G = param$G
  
  #inicialindo a cadeia
  chain = matrix(nrow = r, ncol =  N + 1)
  chain[, 1] = theta_current
  
  for(i in 2:(N + 1)){
    #(i) gerando p_current 
    p_current = mvtnorm::rmvnorm( n = 1, sigma = G(theta_current, param) )
    
    #(ii) Generalizaded Leapfrog
    L = sample(min_L:max_L, 1)
    Prop = glf(eps, L, theta_current, p_current, fixed_p, param)
    
    # (iii) Calculando H's
    h = H(p_current, theta_current, param)
    h = h - H(Prop[ , 2], Prop[ , 1], param)
    
    #(iv) Prop aceitação
    if( runif(1) < min(1, exp(h)) ){
      chain[, i] = Prop[, 1]
      theta_current = matrix(chain[, i], ncol = 1)
      acc = acc + 1
    }else{
      chain[, i] = theta_current
    }
  }
  time.final = Sys.time()
  return(list("chain" = chain, "acc" = acc, "time" = time.final - time.init))
}


mu = 0
sigma = sqrt(2)
T = 3e3
y = rnorm(T, mu, sigma)

N = 1e3
z = rmhmc_R(N, eps = c(0.1, 0.1), min_L = 20, max_L = 30, 
            fixed_p = 5, theta_init = c( -10, 0.5 * log( 10 ) ), 
            param = list(y = y, 
                         mu_0 = 10, s_0 = 3.2, 
                         a = 2, b = 2,
                         logpost = logpost_theta,
                         glogpost = glogpost_theta,
                         G = G_theta_t,
                         d_G = c(d_G_mu, d_G_sigma2),
                         log_J = log_J_theta ),
            #seed = seed
)
z$acc/N
z$time

chain = unlist(z$chain)
# sigma2 = exp(2w)
chain[2, ] = exp( 2 * chain[2, ] )
plot(chain[1, ], chain[2, ], xlab = 'mu', ylab = 'sigma2')

############################### Convergence analysis
# Trace plots
burn = 100
burned = chain[, -c(1:burn)]
plot(burned[1, ], burned[2, ], xlab = 'mu', ylab = 'sigma2')
par(mfrow = c(2, 2))
plot(burned[1, ], type = 'l')
plot(acf(burned[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[2, ], type = 'l')
plot(acf(burned[2, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

# Jumps
lags = 5
jumps = seq(0, N - burn, by = 1 + lags)
burned_lag = burned[, jumps]
par(mfrow = c(2, 2))
plot(burned_lag[1, ], type = 'l')
plot(acf(burned_lag[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[2, ], type = 'l')
plot(acf(burned_lag[2, ], lag.max = 100, plot = FALSE)[1:100])
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
N_new = ncol( burned_lag )
N_eff = coda::effectiveSize(mcmcchain)
IF = N_new / N_eff
IF
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error = round( apply( burned_lag, MARGIN = 1, FUN = sd) / sqrt( N_eff ), 5 )
mc_error

par(mfrow = c(2, 3))
plot(burned_lag[1, ], type = 'l')
plot(acf(burned_lag[1, ], lag.max = 100, plot = FALSE)[1:100])
hist(burned_lag[1, ], breaks = 30, xlab = '', main = '')
plot(burned_lag[2, ], type = 'l')
plot(acf(burned_lag[2, ], lag.max = 100, plot = FALSE)[1:100])
hist(burned_lag[2, ], breaks = 30, xlab = '', main = '')
par(mfrow = c(1, 1))

data = matrix(c(mean(burned_lag[1, ]), 
                quantile(burned_lag[1, ], probs = c(0.025, 0.975)),
                mean(burned_lag[2, ]),
                quantile(burned_lag[2, ], probs = c(0.025, 0.975))),
              nrow = 2, byrow = TRUE)
row.names(data) = c('mu_hat', 'sigma2_hat')
colnames(data) = c('mean', '2.5%', '97,5%')
data

