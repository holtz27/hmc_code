// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//#############################################################################
//#############################################################################
// ################################ Functions #################################
//########################## theta = (mu, phi, sigma) #########################
//######## Transformação: T(theta) = theta'                         ###########  
//######## theta' = (mu, arctanh(phi), log(sigma)) = (mu, w, gama) ############
double logpost_theta(vec theta, List param){
  //#theta = (mu, w, gama)
  //#param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double L = 0.0; 
  double mu = theta[0];
  double phi = tanh( theta[1] ); 
  double sigma = exp( theta[2] );
  
  double  mu_0 = param["mu_0"], s_0 = param["s_0"], 
    a_phi = param["a_phi"], b_phi = param["b_phi"],
    a_s = param["a_s"], b_s = param["b_s"];
  
  vec h = param["h"]; 
  
  //# construindo o vetor z = ( h[2:T] - phi*h[1:(T-1)] - mu*(1 - phi) );
  int T = h.n_elem;
  //vec z = h.subvec(1, T-1) - phi * h.subvec(0, T-2) - mu * (1 - phi);
  
  L += 0.5 * log(1 - phi * phi ) - T * log( sigma ); 
  L -= 0.5 * (1 - phi * phi ) * ( (h[0] - mu)/sigma ) * ( (h[0] - mu)/sigma ); 
  
  
  for( int j = 1 ; j < T ; j++ ){
    L += - (0.5 / (sigma * sigma)) * ( h[j] - mu - phi * (h[j-1] - mu)) * ( h[j] - mu - phi * (h[j-1] - mu));
  }
  
  
  //L -= ( 0.5 * pow(sigma, -2) ) * dot(z, z); 
  
  
  //# priori mu
  L += - 0.5 * ((mu - mu_0)/s_0) * ((mu - mu_0)/s_0);
  //# priori phi
  L += (a_phi - 1) * log(1 + phi) + (b_phi - 1) * log(1 - phi);
  //# priori sigma
  L += - 2 * (a_s + 1) * log(sigma) - b_s / ( sigma * sigma);
  //# jacobiano de T
  L += log( 1 - phi * phi ) + log( sigma );
  
  return L; 
}
vec glogpost_theta(vec theta, List param){
  //#theta = (mu, w, gama)
  //#param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    mu_0 = param["mu_0"], s_0 = param["s_0"], 
                                     a_phi = param["a_phi"], b_phi = param["b_phi"],
                                                                          a_s = param["a_s"], b_s = param["b_s"];
  
  vec h = param["h"], grad = zeros<vec>(3, 1); 
  
  int T = h.n_elem;
  
  // gradiente mu
  grad[0] += (1 - phi * phi ) * (h[0] - mu) / (sigma * sigma); 
  
  for( int j = 1 ; j < T ; j++ ){
    grad[0] += (1 - phi) / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
  
  //priori
  grad[0] += - (mu - mu_0) / (s_0 * s_0);
  
  //gradiente w
  grad[1] += - phi + phi * (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma);
  
  for( int j = 1 ; j < T ; j++ ){
    grad[1] += (1 - phi * phi) / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) ) * (h[j-1] - mu);
  }
  
  //priori
  grad[1] += (a_phi - 1) * (1 - phi);
  grad[1] += - (b_phi - 1) * (1 + phi);
  //jacobiano
  grad[1] += - 2 * phi;
  
  // gradiente gama
  grad[2] += - T + (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma); 
  
  for( int j = 1 ; j < T ; j++ ){
    grad[2] += 1.0 / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
  
  // priori
  grad[2] += - 2 * (a_s + 1) + 2 * b_s / (sigma * sigma);
  // jacobiano
  grad[2] += 1;
  
  return -grad;
}
mat G_theta(vec theta, List param){
  // theta = (mu, w, gama)
  // param = (h, mu_0, s_0, a_phi, b_phi, a_s, b_s)
  double mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]), 
    s_0 = param["s_0"], a_phi = param["a_phi"], b_phi = param["b_phi"], 
                b_s = param["b_s"];
  
  //phi = tanh(theta[1]), sigma = exp(theta[2]),
  
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
  
  //phi = tanh(theta[1]), sigma = exp(theta[2]), 
  
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
  
  //phi = tanh(theta[1]), sigma = exp(theta[2]), 
  
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
//#############################################################################
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
  
  vec y = param["y_T"], h = param["h"]; //, l = param["l"]
  
  //# construindo o vetor z
  //z = sqrt(l * exp(-h)) * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )
  int T = h.n_elem;
  //vec z = sqrt(l%exp(-h))%( y.subvec(1, T)-b0-b1*y.subvec(0, T-1)-b2*exp(h) ); 
  vec z = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp(h);
  vec u = exp( -h ) % z; //vec u = l % exp( -h ) % z;
  
  L = - 0.5 * dot(z, u);
  //# priori b0  
  L += - 0.5 * ((b0 - mu_b0)/s_b0) * (b0 - mu_b0)/s_b0;
  //# priori b1
  L += (a_b1 - 1) * log(1 + b1) + (b_b1 - 1) * log(1 - b1);
  //# priori b2
  L += - 0.5 * ((b2 - mu_b2)/s_b2) * ((b2 - mu_b2)/s_b2);
  //# jacobiano
  L += log(1 - b1 * b1 );
  
  return L;
}
vec glogpost_b(vec b, List param){
  //#b = (b0, delta, b2)
  //#param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //# y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], 
    mu_b0 = param["mu_b0"], s_b0 = param["s_b0"],
                                        a_b1 = param["a_b1"], b_b1 = param["b_b1"],
                                                                          mu_b2 = param["mu_b2"], s_b2 = param["s_b2"];
  
  vec grad = zeros<vec>(3, 1), y = param["y_T"], h = param["h"];
  
  // construindo os vetores u, v e z
  int T = h.n_elem;
  vec z = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp( h );
  vec v = exp( -h );
  vec u = exp( -h ) % y.subvec(0, T-1);
  
  grad[0] = dot(v, z) - (b0 - mu_b0) / (s_b0 * s_b0);
  
  grad[1] = (1 - b1 * b1 ) * ( dot(u, z) ) + (a_b1 - 1) * (1 - b1);
  grad[1] -=  (b_b1 - 1) * (1 + b1); 
  //# jacobiano
  grad[1] -=  2 * b1;
  
  grad[2] = sum( z ) - (b2 - mu_b2) /(s_b2 * s_b2);
  
  return - grad;
}
mat G_b(vec b, List param){
  //b = (b0, delta, b2)
  //param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], s_b0 = param["s_b0"], 
                                                            a_b1 = param["a_b1"], b_b1 = param["b_b1"], s_b2 = param["s_b2"];
  
  vec y = param["y_T"], h = param["h"]; //, l = param["l"]
  int T = h.n_elem;
  
  mat G(3, 3, fill::zeros);
  
  //1a col
  G(0, 0) = sum( exp( -h ) ) + pow(s_b0, -2);
  G(1, 0) = (1 - pow(b1, 2) ) * sum( exp( -h ) % y.subvec(0, T-1) );
  G(2, 0) = T;
  
  //2a col
  G(0, 1) = G(1, 0);
  G(1, 1) = pow(1 - pow(b1, 2), 2) * sum( exp( -h ) % pow(y.subvec(0, T-1), 2) ); 
  G(1, 1) +=  ( 1 - pow(b1, 2) ) * (a_b1 + b_b1 - 2);
  G(2, 1) = (1 - pow(b1, 2) ) * sum( y.subvec(0, T-1) );
  
  //3a col
  G(0, 2) = G(2, 0);
  G(1, 2) = G(2, 1);
  G(2, 2) = sum( exp( h ) ) + pow(s_b2, -2);
  
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
  
  vec y = param["y_T"], h = param["h"]; //, l = param["l"]
  int T = h.n_elem;
  
  mat dG(3, 3, fill::zeros);
  
  //1a col
  dG(1, 0) = -2 * b1 * (1 - pow(b1, 2) ) * sum( exp( -h ) % y.subvec(0, T-1) );
  
  //2a col
  dG(0, 1) = dG(1, 0);
  dG(1, 1) = -2 * b1 * (1 - pow(b1, 2) ) * ( 2 * (1 - pow(b1, 2) ) * sum( exp( -h ) % pow(y.subvec(0, T-1), 2) ) + a_b1 + b_b1 - 2);
  dG(2, 1) = -2 * b1 * (1 - pow(b1, 2) ) * sum( y.subvec(0, T-1) );   
  
  //3a col
  dG(1, 2) = dG(2, 1);
  
  return dG;
}
mat dG_b_b2(vec b, List param){
  mat dG(3, 3, fill::zeros);
  return dG;
}
double log_Jac_b(vec b, List param){
  
  double b1 = tanh( b[1] );
  
  return log(1 - b1 * b1);
}
//########################## h
double logpost_h(vec h, List param){
  //#h = (h1, ..., hT)
  //#param = (y, l, theta, b)
  
  vec y = param["y_T"], theta = param["theta"], b = param["b"]; 
  double L = 0.0, mu = theta[0], phi = tanh( theta[1] ), sigma = exp(theta[2]),
    b0 = b[0], b1 = tanh( b[1] ), b2 = b[2];
  
 
  int T = h.n_elem;
  
  L +=  - 0.5 * sum( h );
 
  for( int j = 1 ; j < T ; j++ ){
    L += - 0.5 / ( sigma * sigma ) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
    L += - 0.5 * exp( -h[j] ) * ( y[j+1] - b0 - b1 * y[j] - b2 * exp( h[j]) ) * ( y[j+1] - b0 - b1 * y[j] - b2 * exp(h[j]) );
  }
  
  L += - 0.5 * exp( -h[0] ) * ( y[1] - b0 - b1 * y[0] - b2 * exp( h[0]) ) * ( y[1] - b0 - b1 * y[0] - b2 * exp(h[0]) );
  L += - 0.5 * (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma);
  
  return L;
}
vec glogpost_h(vec h, List param){
  //#h = (h1, ..., hT)
  //#param = (y, l, theta, b)
  
  vec y = param["y_T"], theta = param["theta"], b = param["b"]; 
  double mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]),
    b0 = b[0], b1 = tanh( b[1] ), b2 = b[2];
  
  //# construindo o vetor s
  int T = h.n_elem;
  vec r  = zeros<vec>(T, 1);
  vec mu_t = y.subvec(1, T) - b0 - b1 * y.subvec(0, T - 1) - b2 * exp(h);
  
  vec s = - 0.5 + 0.5 * exp( -h ) % mu_t % mu_t + b2 * mu_t;
  
  //# construindo o vetor r
  r[0] += ( h[0] - phi * h[1] - mu * (1 - phi) ) / (sigma * sigma);
  r[T-1] = 1 / (sigma * sigma) * (h[T-1] - mu - phi * (h[T-2] - mu) );
  
  vec u = h.subvec(1, T - 2);
  vec v = h.subvec(2, T - 1) + h.subvec(0, T - 3);
  
  r.subvec(1, T-2) += (1 + phi * phi ) * u - phi * v - mu * (1 - phi) * (1 - phi);
  r.subvec(1, T-2) /= sigma * sigma ;
  
  return s - r;
}
//#############################################################################
//#############################################################################
// ################################ Algorithm #################################

//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);
typedef mat    ( * mat_ptr )(vec, List);

// ############################## hmc function
// energy function H
double H_hmc(vec theta, vec p, List param, mat inv_M, num_ptr fun){
  // M = M.i()
  double u, k;
  vec y = param["y_T"];
  
  k = 0.5 * dot(p, inv_M * p);
  u = - fun(theta, param);
  
  return u + k;
}
// generalized leapfrog function
mat lf(double eps, int L, vec theta_current, vec p_current, mat inv_M, 
       List param, vec_ptr fun){
  // M = M.i()
  int T = theta_current.n_elem;
  vec theta = theta_current, p = p_current, y = param["y_T"];
  mat P(T, 2);
  
  // Integration Leapfrog
  for( int k = 0; k < L; k++ ){
    p += 0.5 * eps * fun(theta, param);
    theta += eps * inv_M * p;
    p += 0.5 * eps * fun(theta, param);
  }
  
  P.col(0) = theta;
  P.col(1) = p;
  
  return P;
}
List hmc(int N, double eps, int min_L, int max_L, vec theta_init, 
         List param, num_ptr fun1, vec_ptr fun2){
  
  wall_clock timer;
  timer.tic();
  
  int acc = 0, a = floor(0.1*N), dim = theta_init.n_elem;
  double h, it = N;
  mat chain(dim, N + 1, fill::value(1.7));
  //inicialindo a cadeia
  chain.col(0) = theta_init;
  
  vec mean(dim, fill::zeros), theta_current = theta_init, p_current;
  mat M(dim, dim, fill::eye), inv_M = M.i(), Proposal(dim, 2);
  
  for(int i = 1; i <= N; i++){
    //(i) gerando p_current e u
    p_current = mvnrnd(mean, M);
    
    int L = randi( distr_param( min_L, max_L ) );
    
    //(ii) Leapfrog 
    Proposal = lf(eps, L, theta_current, p_current, inv_M, param, fun2);
    
    // (iii) Calculando H's
    h = H_hmc(theta_current, p_current, param, inv_M, fun1);
    h -= H_hmc(Proposal.col(0), Proposal.col(1), param, inv_M, fun1);
    
    //(iv) Prop aceitação
    if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
      chain.col(i) = Proposal.col(0);
      theta_current = Proposal.col(0);
      acc += 1;
    }else{
      chain.col(i) = theta_current;
    }
    //Progress
    //if( (i % a) == 0 ) cout << "Progresso em " << ceil(100*i/N)<<" %"<< endl;
  }
  
  double time = timer.toc();
  
  return List::create(Named("chain") = chain, Named("acc") = acc, Named("time") = time);
}
// ############################## hmc function one proposal
List hmc_one(double eps, int min_L, int max_L, vec theta_init, 
             List param, num_ptr fun1, vec_ptr fun2){
  
  int acc = 0, dim = theta_init.n_elem; 
  double h; 
  vec mean(dim, fill::zeros), theta_current = theta_init, p_current;
  mat M(dim, dim, fill::eye), inv_M = M.i(), Proposal(dim, 2);
  
  //(i) gerando p_current e u
  p_current = mvnrnd(mean, M);
  
  //(ii) Leapfrog
  int L = randi( distr_param( min_L, max_L ) );
  Proposal = lf(eps, L, theta_current, p_current, inv_M, param, fun2);
  
  // (iii) Calculando H's
  h = H_hmc(theta_current, p_current, param, inv_M, fun1);
  h -= H_hmc(Proposal.col(0), Proposal.col(1), param, inv_M, fun1);
  
  //(iv) Prop aceitação
  if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
    //chain.col(i) = Proposal.col(0);
    theta_current = Proposal.col(0);
    acc += 1;
  }
  
  return List::create(Named("theta_current") = theta_current, Named("acc") = acc);
}
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
// generalized leapfrog function
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
List rmhmc(int N, vec eps, int min_L, int max_L, vec theta_init, List param, int fixed_p, 
           num_ptr fun1, num_ptr log_Jac, vec_ptr fun2, mat_ptr M, mat_ptr v[], int param_id){
  
  //wall_clock timer;
  //timer.tic();
  
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
    if( M(theta_init, param).is_sympd() ){
      p_current  = mvnrnd(mu, M(theta_init, param) );  
    }else{
      cout << "Error in parameter: " << param_id << "!" << endl;
    }
    
    //(ii) Generalizaded Leapfrog
    int L = randi( distr_param( min_L, max_L ) );
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
  
  //double time = timer.toc();
  
  return( List::create(Named("chain") = chain, Named("acc") = acc) );
}
// ############################## rmhmc function one proposal
List rmhmc_one(vec eps, int min_L, int max_L, vec theta_init, List param, int fixed_p, 
               num_ptr fun1, num_ptr log_Jac, vec_ptr fun2, mat_ptr M, mat_ptr v[], int param_id){
  
  // fun1 = logpost
  // fun2 = glogpost
  // M = G
  // v = dG
  
  int acc = 0, r = theta_init.n_elem;
  double h;
  vec theta_current = theta_init, mu(r, fill::zeros), p_current;
  mat Prop(r, 2);
  
  //(i) gerando p_current
  if( M(theta_init, param).is_sympd() ){
    p_current  = mvnrnd(mu, M(theta_init, param) );  
  }else{
    cout << "Error in parameter: " << param_id << "!" << endl;
  }
  
  //(ii) Generalizaded Leapfrog
  int L = randi( distr_param( min_L, max_L ) );
  Prop = glf(eps, L, theta_current, p_current, fixed_p, param, fun2, M, v);
  
  //(iii) Calculando H's
  h = H(theta_current, p_current, param, fun1, log_Jac, M) 
    - H(Prop.col(0), Prop.col(1), param, fun1, log_Jac, M);
  
  //(iv) Prop aceitação
  if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
    theta_current = Prop.col(0);
    acc += 1;
  }
  
  return( List::create(Named("theta_current") = theta_current, Named("acc") = acc) );
}
// ############################## set seed function
void set_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}
// [[Rcpp::export]]
List normal_svm_smn(int N, int T,
                    vec eps_theta, int min_L_theta, int max_L_theta,
                    vec eps_b, int min_L_b, int max_L_b,
                    double eps_h, int min_L_h, int max_L_h,
                    vec init, List param, int fixed_p = 5){
  
  int seed = param["seed"];
  
  if( seed != 0 ) set_seed( seed );
  
  wall_clock timer;
  timer.tic();
  
  mat_ptr v_theta[3] = {&dG_theta_mu, &dG_theta_phi, &dG_theta_sigma},
    v_b[3] = {&dG_b_b0, &dG_b_b1, &dG_b_b2};
  
  vec acc(3, fill::zeros); 
  int z_acc, a = floor(0.1 * N); 
  List z;
  
  // iniciando a cadeia
  mat chain(T + 6, N + 1, fill::zeros); 
  chain.col(0) += init;
  
  for( int it = 1; it <= N; it++ ){
    
    // (i) theta
    z = rmhmc_one(eps_theta, min_L_theta, max_L_theta, 
                  chain.col(it - 1).subvec(0, 2),
                  param, 
                  fixed_p,  
                  &logpost_theta, &log_Jac_theta, &glogpost_theta, &G_theta, v_theta,
                  10);
    vec pivot_1 = z["theta_current"];
    z_acc =  z["acc"];
    acc(0) += z_acc;
    z_acc = 0;
    chain.col( it ).subvec( 0, 2 ) += pivot_1;
    // Atualizando lista param
    param["theta"] = chain.col( it ).subvec( 0, 2 );
    
    // (ii) b
    z = rmhmc_one(eps_b, min_L_b, max_L_b, 
                  chain.col(it - 1).subvec(3, 5), 
                  param, 
                  fixed_p,  
                  &logpost_b, &log_Jac_b, &glogpost_b, &G_b, v_b,
                  20);
    vec pivot_2 = z["theta_current"];
    z_acc =  z["acc"];
    acc(1) += z_acc;
    z_acc = 0;
    chain.col( it ).subvec( 3, 5 ) += pivot_2;
    // Atualizando lista param
    param["b"] = chain.col( it ).subvec( 3, 5 );
    
    // (iii) h
    z = hmc_one(eps_h, min_L_h, max_L_h, 
                chain.col(it - 1).subvec( 6, T + 5 ), 
                param,
                &logpost_h, &glogpost_h);
    
    //mat pivot_3 = z["chain"];
    vec pivot_3 = z["theta_current"];
    z_acc =  z["acc"];
    acc(2) += z_acc;
    z_acc = 0;
    chain.col( it ).subvec( 6, T + 5 ) += pivot_3;
    // Atualizando lista param
    param["h"] = chain.col( it ).subvec( 6,  T + 5 );
    
    //Progress bar
    if( (it % a) == 0 ) cout << "Progresso em " << ceil(100 * it / N)<<" %"<< endl;
  }
  
  double time = timer.toc();
  
  return List::create( Named("chain") = chain, Named("acc") = acc, Named("time") = time );
}
