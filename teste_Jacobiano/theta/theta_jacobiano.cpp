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
  
  vec h = param["h"]; 
  
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
  //L += log( 1 - pow(phi, 2) ) + log(sigma);
  
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
  grad[0] += (1 - phi) * pow(sigma, -2) * sum( z );
  //priori
  grad[0] += - (mu - mu_0) * pow(s_0, -2);
  
  //gradiente w
  grad[1] = - phi + phi * (1 - pow(phi, 2) ) * pow( (h[0] - mu)/sigma, 2 );
  grad[1] += (1 - pow(phi, 2) ) * pow(sigma, -2) * dot(z, s); 
  //priori
  grad[1] += (a_phi - 1) * (1 - phi);
  grad[1] += - (b_phi - 1) * (1 + phi);
  //jacobiano
  //grad[1] += - 2 * phi;
  
  // gradiente gama
  grad[2] = - T + (1 - pow(phi, 2)) * pow( (h[0] - mu)/sigma, 2 ); 
  grad[2] += ( dot(z, z) ) * pow(sigma, -2);   
  // priori
  grad[2] += - 2 * (a_s + 1) + 2 * b_s * pow(sigma, -2);
  // jacobiano
  //grad[2] += 1;
  
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
  G(1, 1) = 2 * pow(phi, 2) + (1 - pow(phi, 2) ) * (T - 1 + a_phi + b_phi - 2);
  G(2, 1) = 2 * phi;
  
  //3a col
  //G(0, 2) = G(2, 0);
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
  
  dG(0, 0) = -2 * (1 - pow(phi, 2) ) * pow(sigma, -2) * (phi + (T - 1) * (1 - phi) ) ;
  dG(1 ,1) = 2 * phi * (1 - pow(phi, 2) ) * (4 - (T-1) - a_phi - b_phi );
  dG(1, 2) = 2 * (1 - pow(phi, 2) );
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
double log_Jac_theta(vec theta, List param){
  
  double phi = tanh( theta[1] ), sigma = exp( theta[2] );
  
  return log(1 - phi * phi) + log(sigma);
}
//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);
typedef mat    ( * mat_ptr )(vec, List);

// ############################## rmhmc function
// energy function H
double H(vec theta, vec p, List param, num_ptr fun, num_ptr log_Jac, mat_ptr M){
  
  double E;
  
  E = - fun(theta, param) - log_Jac(theta, param);
  E += 0.5 * log( det( M(theta, param) ) );
  E += 0.5 * dot(p, M(theta, param).i() * p );
  
  return E;  
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
           int fixed_p, 
           num_ptr fun1, num_ptr log_Jac, vec_ptr fun2, mat_ptr M, mat_ptr v[]){
  
  wall_clock timer;
  timer.tic();
  
  // fun1 = logpost
  // fun2 = glogpost
  // M = G
  // v = dG
  
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
    h = H(theta_current, p_current, param, fun1, log_Jac, M) 
      - H(Prop.col(0), Prop.col(1), param, fun1, log_Jac, M);
    
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
// ############################## rmhmc function one proposal
List rmhmc_one(vec eps, int min_L, int max_L, vec theta_init, List param, 
               int fixed_p, 
               num_ptr fun1, num_ptr log_Jac, vec_ptr fun2, mat_ptr M, mat_ptr v[]){
  
  wall_clock timer;
  timer.tic();
  
  // fun1 = logpost
  // fun2 = glogpost
  // M = G
  // v = dG
  
  int N = 1, i = 1;
  int acc = 0, r = theta_init.n_elem;
  double h;
  vec theta_current = theta_init, mu(r, fill::zeros), p_current;
  mat Prop(r, 2); //chain(r, N + 1)
  
  //(i) gerando p_current
  p_current  = mvnrnd( mu, M(theta_init, param) );
  
  int L = randi( distr_param( min_L, max_L ) );
  
  //(ii) Generalizaded Leapfrog
  Prop = glf(eps, L, theta_current, p_current, fixed_p, param, fun2, M, v);
  
  //(iii) Calculando H's
  h = H(theta_current, p_current, param, fun1, log_Jac, M) 
    - H(Prop.col(0), Prop.col(1), param, fun1, log_Jac, M);
  
  //(iv) Prop aceitação
  if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
    theta_current = Prop.col(0);
    acc += 1;
  }
  
  double time = timer.toc();
  
  return( List::create(Named("theta_current") = theta_current, Named("acc") = acc, Named("time") = time) );
}
// ############################## set seed function
void set_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}
// [[Rcpp::export]]
List rmhmc(int N, vec eps, int min_L, int max_L, vec theta_init, List param, 
           int fixed_p){
  
  //if( is.null( param["seed"] ) ) 
  set_seed( param["seed"] );
  
  // one draw
  /*
  mat_ptr v[3] = {&dG_theta_mu, &dG_theta_phi, &dG_theta_sigma};
  List out = rmhmc_one(eps, min_L, max_L, theta_init, param, fixed_p,  
                       &logpost_theta, &log_Jac_theta, &glogpost_theta, &G_theta, v);
  */
  
  // multiple drws
  mat_ptr v[3] = {&dG_theta_mu, &dG_theta_phi, &dG_theta_sigma};
  List out = Rmhmc(N, eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost_theta, &log_Jac_theta, &glogpost_theta, &G_theta, v);
  
  return out;
  
}
