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
double log_gammapdf(double x, double a, double b){
  return log(pow(b, a)*pow(x, (a-1))*exp(-b*x)/tgamma(a));
}
double logpost(vec theta, List param){
  // Transformation sigma = exp(w)
  // theta = (mu, w)
  // mu = theta[0]
  // w = theta[1]
  double x, mu_0 = param["mu_0"], s_0 = param["s_0"];
  vec y = param["y"], Y = log_normpdf( y, theta[0], exp(theta[1]) );
  
  x = sum( Y );
  x += log_normpdf( theta[0], mu_0, s_0 );
  x += log_gammapdf( exp(2 * theta[1]), param["a"], param["b"] );
  x += 2 * theta[1];
  
  return x;
}
vec glogpost(vec theta, List param){
  
  int T = theta.n_elem;
  double mu = theta[0], w = theta[1], mu_0 = param["mu_0"], s_0 = param["s_0"], a = param["a"], b = param["b"];
  vec y = param["y"], grad(T);
  int N = y.n_elem;
  
  // grad_mu
  grad[0] = sum(y - mu) * exp(-2 * w) - (mu - mu_0) * pow(s_0, -2);
  
  grad[1] = - N + sum( pow( ( y - mu ) , 2) ) * exp( - 2 * w );
  grad[1] += 2 * ( a - b * exp( 2 * w ) );
  
  return grad ;
}
mat G(vec theta, List param){
  vec y = param["y"];
  int T = y.n_elem;
  double g11, g22, s_0 = param["s_0"];
  
  g11 =  T * exp(- 2 * theta[1]) - pow(s_0,-2);
  g22 = 2 * T  + 4 * exp(2 * theta[1]);
  
  mat G = mat{ {g11, 0}, {0, g22} };
  
  return( G );
}
mat dG_mu(vec theta, List param){
  mat M(2, 2, fill::zeros);
  return M;  
}
mat dG_sigma2(vec theta, List param){
  vec y = param["y"];
  int T = y.n_elem;
  double g11, g22, b = param["b"];
  
  g11 = -2 * T * exp(- 2 * theta[1]);
  g22 = 8 * b * exp(2 * theta[1]);
  
  mat G = mat{ {g11, 0}, {0, g22} };
  
  return( G );
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
  vec theta_current = theta_init, mu(2, fill::zeros), p_current;
  mat Prop(r, 2), chain(r, N + 1);
  
  //inicialindo a cadeia
  chain.col(0) = theta_current;
  
  for(int i = 1; i <= N; i++){
    //(i) gerando p_current
    p_current  = mvnrnd(mu, G(theta_init, param));
    
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
  
  mat_ptr v[2] = {&dG_mu, &dG_sigma2};
  List out = Rmhmc(N, eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost, &glogpost, &G, v);
  
  return out;
  
}
