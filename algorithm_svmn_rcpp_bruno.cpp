// [[Rcpp::depends( RcppArmadillo )]]

#include "svmn_rcpp_bruno.h"

// [[Rcpp::export]]
List svmn(int N, 
           int L_theta, double eps_theta,
           int L_b, double eps_b, 
           int L_h, double eps_h,
           vec y_T
           ){
  
  wall_clock timer;
  timer.tic();
  
  int T = y_T.n_elem - 1, a = floor( 0.1 * N );
  
  // iniciando theta
  int acc_theta = 0;
  vec theta_cur = zeros<vec>(3, 1);
  theta_cur[ 0 ] += 0.005;
  theta_cur[ 1 ] += 0.5 * ( log( 1 + 0.98 ) - log( 1 - 0.98 ) );
  theta_cur[ 2 ] += log( sqrt( 0.017 ) );
  
  // iniciando h
  int acc_b = 0;
  vec h_cur = zeros<vec>(T, 1);
  h_cur[ 0 ] += 0.005 + sqrt( 0.03 ) / (1 - 0.95 * 0.95 ) * randn();
  for( int kt = 1 ; kt < T ; kt++ ){
    h_cur[ kt ] += 0.005 + 0.95 * ( h_cur[ kt - 1 ] -0.005 ) + sqrt( 0.03 ) * randn();
  }
  
  // iniciando b
  int acc_h = 0;
  vec b_cur = zeros<vec>(3, 1);
  b_cur[ 0 ] += 0.3;
  b_cur[ 1 ] += 0.5 * ( log( 1 + 0.03 ) - log( 1 - 0.03 ) );
  b_cur[ 2 ] += -0.025;
  
  // iniciando cadeia
  mat chain_theta = zeros<mat>( 3, N + 1 );
  chain_theta.col( 0 ) += theta_cur;
  
  mat chain_b = zeros<mat>( 3, N + 1 );
  chain_b.col( 0 ) += b_cur;
  
  mat chain_h = zeros<mat>( T, N + 1 );
  chain_h.col( 0 ) += h_cur;
  
  // chain builting  
  for(int it = 1 ; it < N + 1 ; it ++){
    
    theta_cur = rmhmc_theta( theta_cur, h_cur, 5, L_theta, eps_theta, T, acc_theta );
    b_cur = rmhmc_b( b_cur, h_cur, 5, L_b, eps_b, T, y_T , acc_b );
    h_cur = hmc_h( h_cur, theta_cur, b_cur, L_h, eps_h, T, y_T, acc_h );
    
    // chain update 
    chain_theta.col( it ) += theta_cur;
    chain_b.col( it ) += b_cur;
    chain_h.col( it ) += h_cur;
    
    //Progress
    if( (it % a) == 0 ) cout << "Progresso em " << ceil( 100 * it / N ) <<" %"<< endl;
  }
  
  List chain = List::create( Named("chain_theta") = chain_theta,
                             Named("chain_b") = chain_b,
                             Named("chain_h") = chain_h );
  
  vec acc = zeros<vec>(3, 1);
  acc[ 0 ] += acc_theta;
  acc[ 1 ] += acc_b;
  acc[ 2 ] += acc_h;
  
  double time = timer.toc();
   
  return List::create( Named("chain") = chain, 
                       Named("acc") = acc, 
                       Named("time") = time 
                       ); 
  
}
