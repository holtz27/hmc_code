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
//########################## h
double logpost_h(vec h, List param){
  //#h = (h1, ..., hT)
  //#param = (y, l, theta, b)
  
  vec y = param["y"], l = param["l"], theta = param["theta"], b = param["b"];
  double L, mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]),
    b0 = b[0], b1 = b[1], b2 = b[2];
  
  //# construindo o vetor z
  int T = h.n_elem;
  vec z(T), u(T - 1);
  
  u = h.subvec(1, T-1) - phi * h.subvec(0, T-2) - mu * (1 - phi) ;
  z = sqrt(l % exp(-h) ) % ( y.subvec(1, T) - b0 - b1*y.subvec(0, T-1) - b2 * exp(h) );
  
  L = sum( log(l) ) - sum( h ) - dot(z, z);
  L += - ( dot(u, u) ) * pow(sigma, -2);
  L += - (1 - pow(phi, 2) ) * pow( (h[0] - mu)/sigma, 2 );
  
  return 0.5 * L;
}
vec glogpost_h(vec h, List param){
  //#h = (h1, ..., hT)
  //#param = (y, l, theta, b)
  
  vec y = param["y"], l = param["l"], theta = param["theta"], b = param["b"];
  double L, mu = theta[0], phi = tanh(theta[1]), sigma = exp(theta[2]),
    b0 = b[0], b1 = b[1], b2 = b[2];
  
  //# construindo o vetor s
  int T = h.n_elem;
  vec s(T), r(T);
  
  s =  0.5 * l % exp(-h) % pow( y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp(h) , 2) 
    + b2 * l % ( y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp(h) ) - 0.5;
  
  //# construindo o vetor r
  r[0] = ( h[0] - phi * h[1] - mu * (1 - phi) ) * pow(sigma, -2);
  r[T-1] = ( h[T-1] - phi * h[T-2] + mu * (1 - phi) ) * pow(sigma, -2);
  
  for(int i = 1; i < T - 1 ; i++){
    r[i] = ( (1 + pow(phi, 2) ) * h[i] - phi * ( h[i+1] + h[i-1]) ) * pow(sigma, -2);
  }
  
  return s - r;
}
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
// [[Rcpp::export]]
