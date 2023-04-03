#include "hmc.h"

//##################################### HMC
// Energy function
double H(vec theta, vec p, List param, mat inv_M, num_ptr fun){
  // M = M.i()
  double E;
    
  E = 0.5 * dot(p, inv_M * p);
  E += - fun(theta, param);
  
  return E;
}

// Leapfrog function
mat lf(double eps, int L, vec theta_current, vec p_current, mat inv_M, List param, vec_ptr fun){
  // M = M.i()
  int T = theta_current.n_elem;
  vec theta = theta_current, p = p_current;
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

// hmc function
List Hmc_one(double eps, int min_L, int max_L, vec theta_init, List param, num_ptr fun1, vec_ptr fun2){
  
  //wall_clock timer;
  //timer.tic();
  
  int acc = 0, dim = theta_init.n_elem; //a = floor(0.1*N)
  double h; 
  
  //inicialindo a cadeia
  
  vec mean(dim, fill::zeros), theta_current = theta_init, p_current;
  mat M(dim, dim, fill::eye), inv_M = M.i(), Proposal(dim, 2);
  
  //(i) gerando p_current e u
  p_current = mvnrnd(mean, M);
  
  //(ii) Leapfrog 
  int L = randi( distr_param( min_L, max_L ) );
  Proposal = lf(eps, L, theta_current, p_current, inv_M, param, fun2);
  
  // (iii) Calculando H's
  h = H(theta_current, p_current, param, inv_M, fun1);
  h -= H(Proposal.col(0), Proposal.col(1), param, inv_M, fun1);
  
  //(iv) Prop aceitação
  if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
    theta_current = Proposal.col(0);
    acc += 1;
  }
  
  //Progress
  //if( (i % a) == 0 ) cout << "Progresso em " << ceil(100*i/N)<<" %"<< endl;
  
  
  //double time = timer.toc();
  
  return List::create(Named("theta_current") = theta_current, Named("acc") = acc);
}

//##################################### RMHMC
// energy function H
double H_rmhmc(vec theta, vec p, List param, num_ptr fun, mat_ptr M){
  
  double E;
  
  E = - fun(theta, param) ;
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

// Generalized Leapfrog function
mat glf(vec eps, int L, vec theta_current, vec p_current, int fixed_p, List param, vec_ptr fun, mat_ptr M, mat_ptr v[]){
  
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

// rmhmc function
List Rmhmc_one(vec eps, int min_L, int max_L, vec theta_init, List param, int fixed_p, num_ptr fun1, vec_ptr fun2, mat_ptr M, mat_ptr v[]){
  
  wall_clock timer;
  timer.tic();
  
  // fun1 = logpost
  // fun2 = glogpost
  // M = G
  // v = dG
  
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
  h = H_rmhmc(theta_current, p_current, param, fun1, M) - H_rmhmc(Prop.col(0), Prop.col(1), param, fun1, M);
    
  //(iv) Prop aceitação
  if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
    theta_current = Prop.col(0);
    acc += 1;
  }
    
  double time = timer.toc();
  
  return( List::create(Named("theta_current") = theta_current, Named("acc") = acc, Named("time") = time) );
}
