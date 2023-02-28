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
###################################################################
############## Transformation sigma = exp(w) ######################
# theta = (mu, w)
logpost = function(theta, param){
  mu = theta[1]
  w = theta[2]
  
  x = sum(dnorm(param$y, mu, exp(w), log = TRUE))
  x = x + dnorm(mu, param$mu_0, param$s_0, log = TRUE)
  x = x + dgamma(exp(2*w), shape = param$a, rate = param$b, log = TRUE)
  x = x + 2*w
  return(x)
}
glogpost = function(theta, param){
  mu = theta[1]
  w = theta[2]
  N = length(y)
  
  grad_mu = sum(param$y - mu)*exp(-2*w)
  grad_mu = grad_mu - (mu - param$mu_0)/(param$s_0)**2
  
  grad_w = -N + sum((param$y - mu)**2)*exp(-2*w)
  grad_w = grad_w + 2*(param$a - param$b*exp(2*w))
  
  return( matrix(c(grad_mu, grad_w), nrow = 2) )
}
###################################################################
############## Gibbs Sampler scheme ###############################
# mu
logpost_mu = function(mu, param){
  y = as.vector(param$y)
  v = as.vector(y - mu)
  x = 0.5*sum(v**2)/param$sigma2 
  x = x + 0.5*(mu - param$mu_0)**2/param$s_0**2
  return(-x)
}
glogpost_mu = function(mu, param){
  v = as.vector(param$y - mu)
  x = sum(v)/param$sigma2 - (mu - param$mu_0)/param$s_0**2
  return(x)
}
# sigma2
logpost_s2 = function(sigma2, param){
  T = length(param$y)
  x = - 0.5*T*log(sigma2) - 0.5*sum((y - param$mu)**2)/sigma2
  x = x + (param$a - 1)*log(sigma2) - sigma2*param$b
  return(x)
}
glogpost_s2 = function(sigma2, param){
  T = length(param$y)
  x = - 0.5*T/sigma2  + 0.5*sum((y - param$mu)**2)/sigma2**2
  x = x + (param$a - 1)/sigma2 - param$b
  return(x)
}

###################################################################
lf = function(eps, L, theta_current, p_current, M, param, fun){
  
  p = matrix(p_current, ncol = 1)
  P = matrix(nrow = nrow(p), ncol = 2)
  theta = matrix(theta_current, ncol = 1)
  eps = matrix(eps, ncol = 1)
  
  for(k in 1:L){
    #p = p + 0.5 * eps * glogpost(theta, param)
    p = p + 0.5 * eps * fun(theta, param)
    theta = theta + eps * (solve(M) %*% p);
    #p = p + 0.5 * eps * glogpost(theta, param)
    p = p + 0.5 * eps * fun(theta, param)
  }
  
  P[, 1] = theta;
  P[, 2] = p;
  
  return(P);
}
H = function(p, M, theta, param, fun){
  
  p = matrix(p, ncol = 1)
  
  k = as.numeric(0.5*(t(p) %*% (solve(M) %*% p)))
  #u = - logpost(theta, param)
  u = - fun(theta, param)
  return(u + k)
}
hmc_R = function(N, eps, L, theta_init, param, fun, gfun, seed){
  set.seed(seed)
  
  acc = 0;
  theta_init = matrix(theta_init, ncol = 1)
  r = nrow(theta_init)
  chain = matrix(nrow = r, ncol =  N + 1)
  chain[ , 1] = theta_init
 
  #inicialindo a cadeia
  theta_current = theta_init
  M = diag(1, r)
  
  for(i in 2:(N+1)){
    #(i) gerando p_current e u
    p_current = matrix(rnorm(r), nrow = r)
    p_prop = matrix(0, nrow = r, ncol = 1)
    
    #(ii) Leapfrog 
    Prop = lf(eps, L, theta_current, p_current, M, param, gfun)
    p_prop = Prop[ , 2]
    
    # (iii) Calculando H's
    h = H(p_current, M, theta_current, param, fun)
    h = h - H(p_prop, M, Prop[ , 1], param, fun)
    
    #(iv) Prop aceitação
    if( runif(1) < min(1, exp(h)) ){
      chain[, i] = Prop[, 1]
      theta_current = matrix(chain[, i], nrow = r)
      acc = acc + 1
    }else{
      chain[, i] = theta_current
    }
  }
  
  return(list("chain" = chain, "acc" = acc))
}
