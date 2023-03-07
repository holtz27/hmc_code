// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//###################################################################
//########################## Modelo
//# prioris: mu ~ N(m_0, s_0), sigma2 ~ Gama(a_0, b_0) 
//# theta = (mu, sigma2)
//###################################################################
//############## Transformation sigma = exp(w) ######################
//# theta = (mu, w)
//###################################################################
//###################################################################
//########################## b = (b0, b1, b2)           #######################
//######## Transformação: T(theta) = b'                 #######################
//######## b' = (b0, arctanh(b1), b2) = (b0, delta, b2) #######################
double logpost_b(vec b, List param){
  //#b = (b0, delta, b2)
  //#param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //# y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], L,
    mu_b0 = param["mu_b0"], s_b0 = param["s_b0"],
                                        a_b1 = param["a_b1"], b_b1 = param["b_b1"],
                                                                          mu_b2 = param["mu_b2"], s_b2 = param["s_b2"];
  
  vec y = param["y"], h = param["h"], l = param["l"];
  
  //# construindo o vetor z
  //z = sqrt(l * exp(-h)) * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )
  int T = h.n_elem;
  vec z = sqrt(l%exp(-h))%( y.subvec(1, T)-b0-b1*y.subvec(0, T-1)-b2*exp(h) ); 
  
  L = - 0.5 * dot(z, z);
  //# priori b0  
  L += - 0.5 * pow( (b0 - mu_b0)/s_b0, 2 );
  //# priori b1
  L += (a_b1 - 1) * log(1 + b1) + (b_b1 - 1) * log(1 - b1);
  //# priori b2
  L += - 0.5 * pow( (b2 - mu_b2)/s_b2, 2 );
  //# jacobiano
  L += log(1 - pow(b1, 2 ) );
  
  return L;
}
vec glogpost_b(vec b, List param){
  //#b = (b0, delta, b2)
  //#param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //# y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], L,
    mu_b0 = param["mu_b0"], s_b0 = param["s_b0"],
                                        a_b1 = param["a_b1"], b_b1 = param["b_b1"],
                                                                          mu_b2 = param["mu_b2"], s_b2 = param["s_b2"];
  
  vec grad(3), y = param["y"], h = param["h"], l = param["l"];
  
  // construindo os vetores u, v e z
  vec v = l % exp(-h); 
  int T = h.n_elem;
  vec z(T), u(T);
  
  u = l % exp(- h ) % y.subvec(0, T-1);
  z = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp( h );
  
  
  
  
  grad[0] = dot(v, z) - (b0 - mu_b0) * pow(s_b0, -2);
  
  grad[1] = (1 - pow(b1, 2) ) * ( dot(u, z) ) + (a_b1 - 1) * (1 - b1);
  grad[1] += (b_b1 - 1) * (1 + b1) - 2 * b1;
  
  grad[2] = dot(l, z) - (b2 - mu_b2) * pow(s_b2, -2);
  
  return grad;
}
mat G_b(vec b, List param){
  //b = (b0, delta, b2)
  //param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], s_b0 = param["s_b0"], 
                                                            a_b1 = param["a_b1"], b_b1 = param["b_b1"], s_b2 = param["s_b2"];
  
  vec y = param["y"], h = param["h"], l = param["l"];
  int T = h.n_elem;
  
  mat G(3, 3, fill::zeros);
  
  //1a col
  G(0, 0) = sum( l % exp( -h ) ) + pow(s_b0, -2);
  G(1, 0) = (1 - pow(b1, 2) ) * sum( l % exp( -h ) % y.subvec(0, T-1) );
  G(2, 0) = sum( l );
  
  //2a col
  G(0, 1) = G(1, 0);
  G(1, 1) = pow(1 - pow(b1, 2), 2) * sum( l % exp( -h ) % pow(y.subvec(0, T-1), 2) ) 
    + ( 1 - pow(b1, 2) ) * (a_b1 + b_b1);
  G(2, 1) = (1 - pow(b1, 2) ) * sum( l % y.subvec(0, T-1) );
  
  //3a col
  G(0, 2) = G(2, 0);
  G(1, 2) = G(2, 1);
  G(2, 2) = sum( l % exp( h ) ) + pow(s_b2, -2);
  
  return G;
}
mat dG_b_b0(vec b, List param){
  mat dG(3, 3, fill::zeros);
  return dG;
}
mat dG_b_b1(vec b, List param){
  //b = (b0, delta, b2)
  //param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], s_b0 = param["s_b0"], 
                                                            a_b1 = param["a_b1"], b_b1 = param["b_b1"], s_b2 = param["s_b2"];
  
  vec y = param["y"], h = param["h"], l = param["l"];
  int T = h.n_elem;
  
  mat dG(3, 3, fill::zeros);
  
  //1a col
  dG(1, 0) = -2 * b1 * (1 - pow(b1, 2) ) * sum( l % exp( -h ) % y.subvec(0, T-1) );
  
  //2a col
  dG(0, 1) = dG(1, 0);
  dG(1, 1) = -2 * b1 * (1 - pow(b1, 2) ) * ( 2 * (1 - pow(b1, 2) ) * sum( l % exp( -h ) % pow(y.subvec(0, T-1), 2) ) - a_b1 - b_b1 );
  dG(2, 1) = -2 * b1 * (1 - pow(b1, 2) ) * sum( l % y.subvec(0, T-1) );   
  
  //3a col
  dG(1, 2) = dG(2, 1);
  
  return dG;
}
mat dG_b_b2(vec b, List param){
  mat dG(3, 3, fill::zeros);
  return dG;
}

//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);
typedef mat    ( * mat_ptr )(vec, List);

// ############################## rmhmc function
// energy function H
double H(vec theta, vec p, List param, num_ptr fun, mat_ptr M){
  
  double x;
  
  x = - fun(theta, param) ;
  x += 0.5 * log( det( M(theta, param) ) );
  x += 0.5 * dot(p, M(theta, param).i() * p );
  
  return x;  
}
// grad energy function H
vec dH_p(vec theta, vec p, List param, mat_ptr M){
  return M(theta, param).i() * p;
}
vec dH_theta(vec theta, vec p, List param, vec_ptr fun, mat_ptr M,
             mat_ptr v[]){
  // param = (d_G = (d_G_theta1, ..., d_G_thetaN), ... )
  
  int D = p.n_elem;
  vec u(D), d; 
  mat i_G = M(theta, param).i(), m;
  
  for(int i = 0; i < D; i++){
    m = i_G * v[i](theta, param);
    u[i] = 0.5 * ( trace( m ) - dot( p, m * i_G * p ) );
  }
  
  //d_H_theta
  d = - fun(theta, param) + u; 
  
  return d ;
}
mat glf(vec eps, int L, vec theta_current, vec p_current, int fixed_p, 
        List param, vec_ptr fun, mat_ptr M, mat_ptr v[]){
  
  int D = p_current.n_elem;
  vec p_n = p_current, theta_n = theta_current, theta_hat = theta_current, p_hat;
  mat Proposal(D, 2); 
  
  for(int t = 0; t < L; t++){
    
    p_hat = p_n;
    
    for(int i = 0; i < fixed_p; i++){
      p_hat = p_n - 0.5 * eps % dH_theta(theta_n, p_hat, param, fun, M, v);
    }
    
    for(int i = 0; i < fixed_p; i++){
      theta_hat = theta_n + 0.5 * eps % ( dH_p(theta_n, p_hat, param, M) + dH_p(theta_hat, p_hat, param, M) );
    }
    
    theta_n = theta_hat;
    p_n = p_hat - 0.5 * eps % dH_theta(theta_n, p_hat, param, fun, M, v);
  }
  
  Proposal.col(0) = theta_n;
  Proposal.col(1) = p_n;
  
  return Proposal;
}

List Rmhmc(int N, vec eps, int min_L, int max_L, vec theta_init, List param, 
           int fixed_p, num_ptr fun1, vec_ptr fun2, mat_ptr M, mat_ptr v[]){
  
  wall_clock timer;
  timer.tic();
  
  // fun1 = logpost
  // fun2 = glogpost
  // M = G
  // v = dG
  
  /*
   if(seed == 0){
   arma_rng::set_seed_random();
   }else{
   arma_rng::set_seed(seed);
   }
   */
  
  int acc = 0, r = theta_init.n_elem;
  double h;
  vec theta_current = theta_init, mu(r, fill::zeros), p_current;
  mat Prop(r, 2), chain(r, N + 1);
  
  //inicialindo a cadeia
  chain.col(0) = theta_current;
  
  for(int i = 1; i <= N; i++){
    //(i) gerando p_current
    p_current  = mvnrnd(mu, M(theta_init, param) );
    
    int L = randi( distr_param( min_L, max_L ) );
    
    //(ii) Generalizaded Leapfrog
    Prop = glf(eps, L, theta_current, p_current, fixed_p, param, fun2, M, v);
    
    //(iii) Calculando H's
    h = H(theta_current, p_current, param, fun1, M) 
      - H(Prop.col(0), Prop.col(1), param, fun1, M);
    
    //(iv) Prop aceitação
    if( R::runif(0, 1) < std::min(1.0, exp(h)) ){
      chain.col(i) = Prop.col(0);
      theta_current = Prop.col(0);
      acc += 1;
    }else{
      chain.col(i) = theta_current;
    }
  }
  
  double time = timer.toc();
  
  return( List::create(Named("chain") = chain, Named("acc") = acc, Named("time") = time) );
}

// [[Rcpp::export]]
List rmhmc(int N, vec eps, int min_L, int max_L, vec theta_init, List param, 
           int fixed_p){
  
  mat_ptr v[3] = {&dG_b_b0, &dG_b_b1, &dG_b_b2};
  List out = Rmhmc(N, eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost_b, &glogpost_b, &G_b, v);
  
  return out;
  
}
######################################################################################################

path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/b_rmhmc.cpp'
Rcpp::sourceCpp(path)

N = 1e4
z = rmhmc(N, eps = c(0.1, 0.1, 0.1), 
            min_L = 10, max_L = 80, 
            theta_init = c(mean(y), atanh(0.8), 2), 
            param = list(y = c(y0, y), 
                         h = h, l = l, 
                         mu_b0 = 0, s_b0 = 10, 
                         a_b1 = 5, b_b1 = 1.5, 
                         mu_b2 = 0, s_b2 = 10), 
            fixed_p = 5)
z$acc/N
z$time

chain = unlist(z$chain)
chain[2, ] = tanh(chain[2, ])

#b0 = 0.001, b1 = -0.5, b2 = -0.2
############################### Convergence analysis
# Trace plots
burn = N/2
burned = chain[, -c(1:burn)]
plot(burned[1, ], burned[2, ], xlab = 'b0', ylab = 'b1')
plot(burned[1, ], burned[3, ], xlab = 'b0', ylab = 'b2')
plot(burned[2, ], burned[3, ], xlab = 'b1', ylab = 'b2')

par(mfrow = c(2, 2))
plot(burned[1, ], type = 'l')
plot(acf(burned[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[2, ], type = 'l')
plot(acf(burned[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[3, ], type = 'l')
plot(acf(burned[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

# Jumps
lags = 1
jumps = seq(0, N - burn, by = 1 + lags)
burned_lag = burned[, jumps]

par(mfrow = c(2, 2))
plot(burned_lag[1, ], type = 'l')
plot(acf(burned_lag[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[2, ], type = 'l')
plot(acf(burned_lag[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned_lag[3, ], type = 'l')
plot(acf(burned_lag[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

plot(burned_lag[1, ], burned_lag[2, ], xlab = 'b0', ylab = 'b1')
plot(burned_lag[1, ], burned_lag[3, ], xlab = 'b0', ylab = 'b2')
plot(burned_lag[2, ], burned_lag[3, ], xlab = 'b1', ylab = 'b2')

par(mfrow = c(1, 3))
hist(burned_lag[1, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[2, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[3, ], breaks = 30, xlab = '', main = '')
par(mfrow = c(1, 1))

data = matrix(c( mean(burned_lag[1, ]), 
                 quantile(burned_lag[1, ], probs = c(0.025, 0.975)),
                 mean(burned_lag[2, ]),
                 quantile(burned_lag[2, ], probs = c(0.025, 0.975)),
                 mean(burned_lag[3, ]),
                 quantile(burned_lag[3, ], probs = c(0.025, 0.975)) ), 
              nrow = 3, byrow = TRUE)
row.names(data) = c('b0_hat', 'b1_hat', 'b2_hat')
colnames(data) = c('mean', '5%', '97,5%')
data









