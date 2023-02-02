########################## theta = (mu, phi, sigma) #########################
######## Transformação: T(theta) = theta'                         ############  
######## theta' = (mu, arctanh(phi), log(sigma)) = (mu, w, gama) ############
logpost_theta = function(theta, param){
  #theta = (mu, w, gama)
  #param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  h = as.vector(param$h) #obs. transformar em vetor
  
  # construindo o vetor z
  T = length(h)
  z = ( h[2:T] - phi*h[1:(T-1)] - mu*(1 - phi) )
  
  #L = 0.5 * log(1 - phi**2) - 0.5 * (1 - phi**2)*( (h[1] - mu)/sigma )**2
  #L = L - (0.5/sigma**2) * t(z) %*% z - T * log(sigma)
  # priori mu
  #L = L - 0.5 * ( (mu - param$mu_0)/param$s_0 )**2
  # priori phi
  #L = L + (param$a_phi - 1) * log(1 + phi) + (param$b_phi - 1) * log(1 - phi)
  # priori sigma
  #L = L - 2 * (param$a_s + 1) * log(sigma) - param$b_s/sigma**2
  # jacobiano de T
  #L = L + log(1 - phi**2) + log(sigma)
  
  
  L = ( 
    0.5 * log(1 - phi**2) - 0.5 * (1 - phi**2)*( (h[1] - mu)/sigma )**2
    - (0.5/sigma**2) * t(z) %*% z - T * log(sigma)
    # priori mu
    - 0.5 * ( (mu - param$mu_0)/param$s_0 )**2
    # priori phi
    + (param$a_phi - 1) * log(1 + phi) + (param$b_phi - 1) * log(1 - phi)
    # priori sigma
    - 2 * (param$a_s + 1) * log(sigma) - param$b_s/sigma**2
    # jacobiano de T
    + log(1 - phi**2) + log(sigma)
  )
  
  
  return( as.numeric(L) )
}
glogpost_theta = function(theta, param){
  #theta = (mu, w, gama)
  #param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  
  h = as.vector(param$h) #obs. transformar em vetor
  
  # construindo o vetor z
  T = length(h)
  z = ( h[2:T] - phi * h[1:(T-1)] - mu * (1 - phi) )
  s = h[1:(T-1)] - mu
  
  grad = matrix(nrow = 3, ncol = 1)
  
  # gradiente mu
  grad[1] = (1 - phi**2) * (h[1] - mu) / sigma**2
  grad[1] = grad[1] + (1 - phi)/ sigma * sum( z )
  # priori
  grad[1] = grad[1] - (mu - param$mu_0)/param$s_0**2
  
  # gradiente w
  grad[2] = - phi + phi * (1 - phi**2) * ( (h[1] - mu)/sigma )**2
  grad[2] = grad[2] + (1 - phi**2)/sigma**2 * ( z %*% s ) 
  # priori
  grad[2] = grad[2] + (param$a_phi - 1) * (1 - phi)
  grad[2] = grad[2] - (param$b_phi - 1) * (1 + phi)
  # jacobiano
  grad[2] = grad[2] - 2 * phi
    
  # gradiente gama
  grad[3] = - T + (1 - phi**2) * ( (h[1] - mu)/sigma )**2 + ( t(z) %*% z )/sigma
  # priori
  grad[3] = grad[3] - 2 * (param$a_s + 1) + 2 * param$b_s/sigma**2
  # jacobiano
  grad[3] = grad[3] + 1

  return( matrix(grad, ncol = 1) )  
}
G_theta = function(theta, param){
  #theta = (mu, w, gama)
  #param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  T = length(param$h)
  
  # 1acol
  m11 = ( (1-phi**2) + (T-1) * (1 - phi)**2 )/sigma**2 + 1/param$s_0**2
  m21 = m31 = 0
  
  #2a col
  m12 = 0
  m22 = 2*phi**2 + (1 - phi**2) * (T - 1 + param$a_phi + param$b_phi)
  m32 = 2 * phi
  
  #3a col
  m13 = m31
  m23 = m32
  m33 = 2 * T + 4 * param$b_s/sigma**2
  
  return( matrix(c( m11, m21, m31,
                    m21, m22, m23,
                    m13, m23, m33), nrow = 3, ncol = 3) )
}
dG_theta_mu = function(theta, param) return( matrix(0, nrow = 3, ncol = 3)  )
dG_theta_phi = function(theta, param){
  #theta = (mu, w, gama)
  #param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  T = length(param$h)
  
  m11 = -2 * (1 - phi**2) * (phi + (T-1) * (1 - phi) )/sigma**2
  m22 = 2 * phi * (1 - phi**2) * (2 - (T-1) - param$a_phi - param$b_phi )
  m23 = 2 * (1 - phi)
  
  return( matrix(c(m11, 0, 0,
                   0, m22, m23,
                   0, m23, 0), nrow = 3, ncol = 3) )
}
dG_theta_sigma = function(theta, param){
  #theta = (mu, w, gama)
  #param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  T = length(param$h)
  
  m11 = - 2 * (1 - phi**2 + (T-1) * (1 - phi)**2 )/sigma**2
  m33 = -8 * param$b_s/sigma**2
  
  return( matrix(c(m11, 0, 0,
                   0, 0, 0,
                   0, 0, m33), nrow = 3, ncol = 3) )
}

########################## b = (b0, b1, b2)           #######################
######## Transformação: T(theta) = b'                 #######################
######## b' = (b0, arctanh(b1), b2) = (b0, delta, b2) #######################
logpost_b = function(b, param){
  #b = (b0, delta, b2)
  #param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  # y = (y0. y1, ..., yT)
  
  b0 = b[1]
  b1 = tanh(b[2])
  b2 = b[3]
  y = as.vector(param$y)
  h = as.vector(param$h)
  l = as.vector(param$l)
  
  # construindo o vetor z
  T = length(h)
  z = sqrt(l * exp(-h)) * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )
  
  L = - 0.5 * t(z) %*% z
  # priori b0  
  L = L - 0.5 * ( (b0 - param$mu_b0)/param$s_b0 )**2
  # priori b1
  L = L + (param$a_b1 - 1) * log(1 + b1) + (param$b_b1 - 1) * log(1 - b1)
  # priori b2
  L = L - 0.5 * ( (b2 - param$mu_b2)/param$s_b2 )**2
  # jacobiano
  L = L + log(1 - b1**2)
  
  return( as.numeric(L) )
}
glogpost_b = function(b, param){
  #b = (b0, delta, b2)
  #param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  # y = (y0. y1, ..., yT)
  
  b0 = b[1]
  b1 = tanh(b[2])
  b2 = b[3]
  y = as.vector(param$y)
  h = as.vector(param$h)
  l = as.vector(param$l)
  
  # construindo os vetores u, v e z
  T = length(h)
  v = l * exp(-h)
  u = l * exp(-h) * y[1:T]
  z = ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )
  
  grad = matrix(nrow = 3, ncol = 1)
  
  grad[1] = ( t(v) %*% z ) - (b0 - param$mu_b0)/param$s_b0**2
  
  grad[2] = (1 - b1**2) * ( t(u) %*% z ) + (param$a_b1 - 1) * (1 - b1)
  grad[2] = grad[2] + (param$b_b1 - 1) * (1 + b1) - 2 * b1
  
  grad[3] = t(l) %*% z - (b2 - param$mu_b2)/param$s_b2**2
  
  return( matrix(grad, ncol = 1) )  
}
G_b = function(b, param){
  #b = (b0, delta, b2)
  #param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  # y = (y0. y1, ..., yT)
  
  b0 = b[1]
  b1 = tanh(b[2])
  b2 = b[3]
  y = as.vector(param$y)
  h = as.vector(param$h)
  l = as.vector(param$l)
  T = length(h)
  
  #1a col
  m11 = sum( l * exp( -h ) ) + 1/param$s_b0**2
  m21 = (1 - b1**2) * sum( l * exp( -h ) * y[1:T] )
  m31 = sum( l )
  
  #2a col
  m12 = m21
  m22 = (1 - b1**2)**2 * sum( l * exp( -h ) * y[1:T]**2 ) + (1 - b1**2) * (param$a_b1 + param$b_b1)
  m32 = (1 - b1**2) * sum( l * y[1:T] )
  
  #3a col
  m13 = m31
  m23 = m32
  m33 = sum( l * exp( h ) ) + 1/param$s_b2**2
  
  return( matrix(c(m11, m21, m31,
                   m12, m22, m32,
                   m13, m23, m33), nrow = 3, ncol = 3) )
  
}
dG_b_b0 = function(b, param) return( matrix(0, nrow = 3, ncol = 3)  )
dG_b_b1 = function(b, param){
  #b = (b0, delta, b2)
  #param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  # y = (y0. y1, ..., yT)
  
  b0 = b[1]
  b1 = tanh(b[2])
  b2 = b[3]
  y = as.vector(param$y)
  h = as.vector(param$h)
  l = as.vector(param$l)
  T = length(h)
  
  #1a col
  m11 = m31 = 0
  m21 = -2 * b1 * (1 - b1**2) * sum( l * exp( -h ) * y[1:T] )
  
  #2a col
  m12 = m21
  m22 = -2 * b1 * (1 - b1**2) * ( 2 * (1 - b1**2) * sum( l * exp( -h ) * y[1:T]**2 ) - param$a_b1 - param$b_b1 )
  m32 = -2 * b1 * (1 - b1**2) * sum( l * y[1:T] )
  
  #3a col
  m13 = m33 = 0
  m23 = m32
  
  return( matrix(c(m11, m21, m31,
                   m12, m22, m32,
                   m13, m23, m33), nrow = 3, ncol = 3) )
  
}
dG_b_b2 = function(b, param) return( matrix(0, nrow = 3, ncol = 3)  )
########################## h
logpost_h = function(h, param){
  #h = (h1, ..., hT)
  #param = (y, l, theta, b)
  
  y = as.vector(param$y)
  l = as.vector(param$l)
  
  theta = as.vector(param$theta)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  
  b = as.vector(param$b)
  b0 = b[1]
  b1 = b[2]
  b2 = b[3]
  
  # construindo o vetor z
  T = length(h)
  z = sqrt(l * exp(-h)) * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )
  u = ( h[2:T] - phi*h[1:(T-1)] - mu*(1 - phi) )
  
  L = sum( log(l) ) - sum( h ) - t(z) %*% z
  L = L - ( t(u) %*% u )/sigma**2
  L = L - (1 - phi**2) * ( (h[1] - mu)/sigma )**2
  
  return( as.numeric(0.5 * L) )
}
glogpost_h = function(h, param){
  #h = (h1, ..., hT)
  #param = (y, l, theta, b)
  
  y = as.vector(param$y)
  l = as.vector(param$l)
  
  theta = as.vector(param$theta)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  
  b = as.vector(param$b)
  b0 = b[1]
  b1 = b[2]
  b2 = b[3]
  
  # construindo o vetor s
  T = length(h)
  s =  0.5 * l * exp(-h) * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )**2
  s = s + b2 * l * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) ) - 0.5
  
  # construindo o vetor r
  delta_1 = ( h[1] - phi * h[2] - mu * (1 - phi) )/sigma**2
  w = ( (1 + phi**2) * h[2:(T-1)] - phi * ( h[3:T] + h[1:(T-2)]) )/sigma**2
  delta_T = ( h[T] - phi * h[T-1] + mu * (1 - phi) )/sigma**2
  
  r = c(delta_1, w, delta_T)
  
  return( matrix(s - r, ncol = 1) )
}
G_h = function(h, param){
  #h = (h1, ..., hT)
  #param = (y, l, theta, b)
  
  y = as.vector(param$y)
  l = as.vector(param$l)
  
  theta = as.vector(param$theta)
  mu = theta[1]
  phi = tanh(theta[2])
  sigma = exp(theta[3])
  
  b = as.vector(param$b)
  b0 = b[1]
  b1 = b[2]
  b2 = b[3]
  
  T = length(h)
  
  G = b2**2 * diag( 0.5 + 1/sigma**2 + l * exp(h) )
  G[abs(row(G) - col(G)) == 1] = ( - phi/sigma**2 )
  
  return( G )
}
# dG_h_ht = diag(0, ..., 0, b2² * l[t] * exp( h[t] ), 0, ..., 0)
#############################################################################
########################## l | v ~ Gama(v/2, v/2)

########### l_t | y, b, h, v ~ Gama( a_t, b_t ), 1 < t < T
# a_t = (v + 1)/2
# b_t =  ( exp(-h_t)(y_t - b0 - b1y_{t-1} - b2exp{h_t})**2 + v )/2
l_gibbs = function(param){
  #param = (y, b, v, h)
  y = as.vector(param$y)
  b = as.vector(param$b)
  v = param$v
  h = as.vector(param$h)
  T = length(h)
  u = exp(-h) * ( y[2:(T+1)] - b[1] - b[2] * y[1:T] - b[3] * exp(h) )**2
  
  l = rgamma( T, shape = 0.5 * (v + 1), rate = 0.5 * (u + v) )
  
  return( as.vector(l) )
}
########################## v.1
######## Transformação: T(v) = e                         ############
######## e = log(v)
logpost_v.1 = function(e, param){
  #param = (l, a_v, b_v)
  v = exp(e)
  l = as.vector(param$l)
  T = length(l)
  
  L = 0.5 * T * v * log( 0.5 * v  ) - T * log( gamma( 0.5 * v )  )
  L = L + 0.5 * v * sum( log(l) - l )
  # priori
  L = L + param$a_v * log(v) - param$b_v * v
  
  return( as.numeric(L) )
}
glogpost_v.1 = function(e, param){
  #param = (l, a_v, b_v)
  v = exp(e)
  l = as.vector(param$l)
  T = length(l)
  
  grad = 0.5 * T * v * log(0.5 * v) - 0.5 * T * v * digamma( 0.5 * v)
  grad = grad + 0.5 * v * sum( log(l) - l ) 
  grad = grad + param$a_v - param$b_v * v + 0.5 * T * v
  
  return( as.numeric(grad) )
}
G_v.1 = function(e, param){
  #param = (l, a_v, b_v)
  v = exp(e)
  l = as.vector(param$l)
  T = length(l)
  
  G = T * v**2 * psigamma(0.5 * v, 1) /4 - v * ( param$b_v - 0.5 * T )
  
  return( G )
}
dG_v.1 = function(e, param){
  #param = (l, a_v, b_v)
  v = exp(e)
  l = as.vector(param$l)
  T = length(l)
  
  dG = 0.5 * T * v**2 * psigamma(0.5 * v, 1) + T * v**3 * psigamma(0.5 * v, 2)/8 - 0.5 * T - param$b_v
  
  return( dG )
}
#############################################################################
########################## l | v ~ Beta(v, 1)

# l_t | y, b, h, v ~ Gama_truncada_{0 < l_t < 1}( a_t, b_t ), 1 < t < T
# a_t = v + 1/2
# b_t =  ( exp(-h_t)(y_t - b0 - b1 * y_{t-1} - b2 * exp{h_t})**2 )/2

########################## v.2
# v.2 | l_t ~ Gama_truncada_{1 < v}( a_t, b_t )
# a_t = T + a_v
# b_t = b_v - sum( log(l)  )

#############################################################################
########################## l | v ~ IGama(v/2, v/2)

########### l_t | y, b, h, v ~ GIGauss( a_t, b_t, p_t ), 1 < t < T
# a_t = (1 - v)/2
# b_t =  ( exp(-h_t)(y_t - b0 - b1 * y_{t-1} - b2 * exp{h_t})**2 )/2
# p_t = v

########################## v.3
######## Transformação: T(v) = e                         ############
######## e = log(v)
logpost_v.3 = function(e, param){
  #param = (l, a_v, b_v)
  v = exp(e)
  l = as.vector(param$l)
  T = length(T)
  
  L = 0.5 * T * v * log( 0.5 * v ) - T * log( gamma(0.5 * v) )
  L = L - 0.5 * v * sum(log(l) - 1/l)
  # priori
  L = L + param$a_v * log(v) - param$b_v * v
  
  return( as.numeric(L) )
}
glogpost_v.3 = function(e, param){
  #param = (l, a_v, b_v)
  v = exp(e)
  l = as.vector(param$l)
  T = length(T)
  
  grad = 0.5 * T * v * log( 0.5 * v ) + 0.5 * T * v
  grad = grad - 0.5 * T * v * digamma( 0.5 * v ) - 0.5 * v * sum(log(l) - 1/l)
  grad = grad + param$a_v - param$b_v * v
  
  return( as.numeric(grad) )
}
