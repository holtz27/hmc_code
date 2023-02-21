// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//########################## theta = (mu, phi, sigma) #########################
//######## Transformação: T(theta) = theta'                         ###########  
//######## theta' = (mu, arctanh(phi), log(sigma)) = (mu, w, gama) ############
double logpost_theta(vec theta, List param){
  //#theta = (mu, w, gama)
  //#param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double L, mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    mu_0 = param["mu_0"], s_0 = param["s_0"], 
                                     a_phi = param["a_phi"], b_phi = param["b_phi"],
                                                                          a_s = param["a_s"], b_s = param["b_s"];
  
  vec h = param["h"]; //obs. transformar em vetor
  
  //# construindo o vetor z = ( h[2:T] - phi*h[1:(T-1)] - mu*(1 - phi) );
  int T = h.n_elem;
  vec z(T - 1);
  
  z = h.subvec(1, T-1) - phi * h.subvec(0, T-2) - mu * (1 - phi);
  
  L = 0.5 * log(1 - pow(phi, 2) ); 
  L += - 0.5 * (1 - pow(phi, 2) ) * pow( (h[0] - mu)/sigma , 2); 
  L += - ( 0.5 * pow(sigma, -2) ) * dot(z, z) - T * log(sigma);
  //# priori mu
  L += - 0.5 * pow((mu - mu_0)/s_0, 2);
  //# priori phi
  L += (a_phi - 1) * log(1 + phi) + (b_phi - 1) * log(1 - phi);
  //# priori sigma
  L += - 2 * (a_s + 1) * log(sigma) - b_s * pow(sigma, -2);
  //# jacobiano de T
  L += log( 1 - pow(phi, 2) ) + log(sigma);
  
  return L; 
}
vec glogpost_theta(vec theta, List param){
  //#theta = (mu, w, gama)
  //#param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double L, mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    mu_0 = param["mu_0"], s_0 = param["s_0"], 
                                     a_phi = param["a_phi"], b_phi = param["b_phi"],
                                                                          a_s = param["a_s"], b_s = param["b_s"];
  
  vec h = param["h"], grad(3); //obs. transformar em vetor
  
  // construindo o vetor z = ( h[2:T] - phi*h[1:(T-1)] - mu*(1 - phi) )
  // e o vetor s = h[1:(T-1)] - mu
  int T = h.n_elem;
  vec z(T - 1), s(T - 1);
  
  z = h.subvec(1, T-1) - phi * h.subvec(0, T-2) - mu * (1 - phi);
  s = h.subvec(0, T-2) - mu;
  
  // gradiente mu
  grad[0] = (1 - pow(phi, 2) ) * (h[0] - mu) * pow(sigma, -2); 
  grad[0] += (1 - phi)/sigma * sum( z );
  //priori
  grad[0] += - (mu - mu_0) * pow(s_0, -2);
  
  //gradiente w
  grad[1] = - phi + phi * (1 - pow(phi, 2) ) * pow( (h[0] - mu)/sigma, 2 );
  grad[1] += (1 - pow(phi, 2) ) * pow(sigma, -2) * dot(z, s); 
  //priori
  grad[1] += (a_phi - 1) * (1 - phi);
  grad[1] += - (b_phi - 1) * (1 + phi);
  //jacobiano
  grad[1] += - 2 * phi;
  
  // gradiente gama
  grad[2] = - T + (1 - pow(phi, 2)) * pow( (h[0] - mu)/sigma, 2 ); 
  grad[2] += ( dot(z, z) )/sigma;   
  // priori
  grad[2] += - 2 * (a_s + 1) + 2 * b_s * pow(sigma, -2);
  // jacobiano
  grad[2] += 1;
  
  return grad;
}
mat G_theta(vec theta, List param){
  // theta = (mu, w, gama)
  // param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    s_0 = param["s_0"], a_phi = param["a_phi"], b_phi = param["b_phi"], 
                b_s = param["b_s"];
  
  vec h = param["h"];
  int T = h.n_elem;
  mat G(3, 3, fill::zeros);
  
  //1acol
  G(0, 0) = ( (1 - pow(phi, 2) ) 
                + (T - 1) * pow(1 - phi, 2) ) * pow(sigma, -2) + pow(s_0, -2);
  //G(1, 0) = G(2, 0) = 0;
  
  //#2a col
  //G(0, 1) = 0;
  G(1, 1) = 2 * pow(phi, 2) + (1 - pow(phi, 2) ) * (T - 1 + a_phi + b_phi);
  G(2, 1) = 2 * phi;
  
  //3a col
  G(0, 2) = G(2, 0);
  G(1, 2) = G(2, 1);
  G(2, 2) = 2 * T + 4 * b_s * pow(sigma, -2);
  
  return G;
}
mat dG_theta_mu(vec theta, List param){
  mat dG(3, 3, fill::zeros);
  return dG;
}
mat dG_theta_phi(vec theta, List param){
  // theta = (mu, w, gama)
  // param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    s_0 = param["s_0"], a_phi = param["a_phi"], b_phi = param["b_phi"], 
                b_s = param["b_s"];
  
  vec h = param["h"];
  int T = h.n_elem;
  mat dG(3, 3, fill::zeros);
  
  dG(0, 0) = -2 * (1 - pow(phi, 2) ) * (phi + (T - 1) * (1 - phi) ) * pow(sigma, -2);
  dG(1 ,1) = 2 * phi * (1 - pow(phi, 2) ) * (2 - (T-1) - a_phi - b_phi );
  dG(1, 2) = 2 * (1 - phi);
  dG(2, 1) = dG(1, 2);
  
  return dG;
}
mat dG_theta_sigma(vec theta, List param){
  // theta = (mu, w, gama)
  // param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    b_s = param["b_s"];
  
  vec h = param["h"];
  int T = h.n_elem;
  mat dG(3, 3, fill::zeros);
  
  dG(0, 0) = -2 * (1 - pow(phi, 2) + (T-1) * pow(1 - phi, 2) ) * pow(sigma, -2);
  dG(2, 2) = -8 * b_s * pow(sigma, -2);
  
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
  
  mat_ptr v[3] = {&dG_theta_mu, &dG_theta_phi, &dG_theta_sigma};
  
  List out = Rmhmc(N, eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost_theta, &glogpost_theta, &G_theta, v);
  
  return out;
  
}


###############################################################################################
path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/Parâmetros/theta/theta_rmhmc.cpp'
Rcpp::sourceCpp(path)

#theta = (mu, w, gama)
#param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)

N = 1e4
z = rmhmc(N, eps = c(0.015, 0.015, 0.015), 
          min_L = 50, max_L = 50, 
          theta_init = c( 1, atanh( 0.86 ), log( 10 ) ), 
          param = list(y = c(y0, y), 
                       h = h, 
                       l = l, 
                       mu_0 = 0, s_0 = 3.2, 
                       a_phi = 20, b_phi = 1.5, 
                       a_s = 2.5, b_s = 0.025), 
          fixed_p = 5)
z$acc/N
z$time

chain = unlist(z$chain)
chain[2, ] = tanh(chain[2, ])
chain[3, ] = exp(chain[3, ])
#mu = -1
#phi = 0.95
#sigma = 0.25
############################### Convergence analysis
# Trace plots
burn = N/2 + 500
burned = chain[, -c(1:burn)]

par(mfrow = c(2, 2))
plot(burned[1, ], type = 'l')
plot(acf(burned[1, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[2, ], type = 'l')
plot(acf(burned[2, ], lag.max = 100, plot = FALSE)[1:100])
plot(burned[3, ], type = 'l')
plot(acf(burned[3, ], lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))
# Jumps
lags = 10
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

plot(burned[1, ], burned[2, ], xlab = 'mu', ylab = 'phi')
plot(burned[1, ], burned[3, ], xlab = 'mu', ylab = 'sigma')
plot(burned[2, ], burned[3, ], xlab = 'phi', ylab = 'sigma')
###########################################################################
N_new = length( burned )
mcmcchain = coda::as.mcmc( t(burned) )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff = coda::effectiveSize(mcmcchain)
IF = N_new / N_eff
IF
# |G| > 1.96 evidencia não convergencia
coda::geweke.diag(mcmcchain)

plot(burned_lag[1, ], burned_lag[2, ], xlab = 'b0', ylab = 'b1')
plot(burned_lag[1, ], burned_lag[3, ], xlab = 'b0', ylab = 'b2')
plot(burned_lag[2, ], burned_lag[3, ], xlab = 'b1', ylab = 'b2')

par(mfrow = c(1, 3))
hist(burned_lag[1, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[2, ], breaks = 30, xlab = '', main = '')
hist(burned_lag[3, ], breaks = 30, xlab = '', main = '')
par(mfrow = c(1, 1))

data = matrix(c( mean(burned_lag[1, ]), 
                 quantile(burned_lag[1, ], probs = c(0.05, 0.975)),
                 mean(burned_lag[2, ]),
                 quantile(burned_lag[2, ], probs = c(0.05, 0.975)),
                 mean(burned_lag[3, ]),
                 quantile(burned_lag[3, ], probs = c(0.05, 0.975)) ), 
              nrow = 3, byrow = TRUE)
row.names(data) = c('b0_hat', 'b1_hat', 'b2_hat')
colnames(data) = c('mean', '5%', '97,5%')
data





















