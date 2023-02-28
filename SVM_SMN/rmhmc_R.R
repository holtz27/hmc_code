########################## Modelo
# prioris: mu ~ N(m_0, s_0), sigma2 ~ Gama(a_0, b_0) 
# theta = (mu, sigma2)
logpost = function(theta, param){
  x = sum(dnorm(param$y, theta[1], sqrt(theta[2]), log = TRUE))
  x = x + dnorm(theta[1], param$mu_0, param$s_0, log = TRUE)
  x = x + dgamma(theta[2], shape = param$a, rate = param$b, log = TRUE)
  return(x)
}
glogpost = function(theta, param){
  N = length(param$y)
  
  grad_mu = sum(param$y-theta[1])/theta[2]-(theta[1]-param$mu_0)/param$s_0**2
  
  grad_sigma2 = 0.5*sum((param$y - theta[1])**2)/theta[2]**2 - param$b
  grad_sigma2 = grad_sigma2 - 0.5*N/theta[2] + (param$a-1)/theta[2]
  
  return( matrix(c(grad_mu, grad_sigma2), nrow = 2) )
}
G = function(theta, param){
  T = length(param$y)
  g11 = T/theta[2] - 1/param$s_0**2
  g22 = T/theta[2]**2 - (param$a - 1)/theta[2]**2
  return(diag(c(g11, g22), 2, 2))
}
d_G_mu = function(theta, param) return(matrix(0, 2, 2))
d_G_sigma2 = function(theta, param){
  T = length(param$y)
  g11 = -T/theta[2]**2
  g22 = -2*( 0.5*T - param$a + 1 )/theta[2]**3
  return(diag(c(g11, g22), 2, 2))
}

########################## Algoritmo

H = function(p, theta, param){
  p = matrix(p, ncol = 1)
  x = - logpost(theta, param) 
  x = x + 0.5 * log( det( G(theta, param) ) )  
  x = x + 0.5 *  t(p) %*% solve( G(theta, param) ) %*% p 
 
 return(x)
}
d_H_theta = function(p, theta, param){
  # param = (d_G = (d_G_theta1, ..., d_G_thetaN), ... )
  p = matrix(p, ncol = 1)
  D = length(p)
  d = matrix(nrow = D, ncol = 1)
  G.i = solve(G(theta, param))
  
  u = matrix(ncol = 1)
  
  for(i in 1:D){
    M = G.i %*% param$d_G[[i]](theta, param)
    u[i] = 0.5 * sum(diag( M ))  
    u[i] = - 0.5 * t(p) %*% M %*% d_H_p(p, theta, param) 
  }
  
  u = matrix(u, ncol = 1)
  
  #d_H_theta
  d = - glogpost(theta, param) + u 
  
  return( d )
}
d_H_p = function(p, theta, param){
  p = matrix(p, ncol = 1)
  return( solve( G(theta, param) ) %*% p )
}
glf = function(eps, L, theta_current, p_current, fixed_p, param){
  
  p_n = matrix(p_current, ncol = 1)
  theta_n = theta_hat = matrix(theta_current, ncol = 1)
  Proposal = matrix(nrow = nrow(p_n), ncol = 2)
  eps = matrix(eps, ncol = 1)
  
  for(t in 1:L){
    
    p_hat = p_n
    
    for(i in 1:fixed_p) p_hat = p_n - 0.5 * eps * d_H_theta(p_hat, theta_n, param)
    
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
rmhmc_R = function(N, eps, min_L, max_L, theta_init, param, fixed_p, seed = NULL){
  
  if ( !is.null(seed) ) set.seed(seed)
  time.init = Sys.time()
  acc = 0;
  theta_current = matrix(theta_init, ncol = 1)
  r = nrow(theta_current)
  
  #inicialindo a cadeia
  chain = matrix(nrow = r, ncol =  N + 1)
  chain[, 1] = theta_current
  
  for(i in 2:(N + 1)){
    #(i) gerando p_current 
    p_current = mvtnorm::rmvnorm(n = 1, sigma = G(theta_current, param))
    #p_prop = matrix(0, nrow = r, ncol = 1)
    
    L = sample(min_L:max_L, 1)
    
    #(ii) Generalizaded Leapfrog
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
