// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
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
  
  vec y = param["y_T"], h = param["h"], l = param["l"];
  
  //# construindo o vetor z
  //z = sqrt(l * exp(-h)) * ( y[2:(T+1)] - b0 - b1 * y[1:T] - b2 * exp(h) )
  int T = h.n_elem;
  //vec z = sqrt(l%exp(-h))%( y.subvec(1, T)-b0-b1*y.subvec(0, T-1)-b2*exp(h) ); 
  vec z = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp(h);
  vec u = l % exp( -h ) % z;
  
  L = - 0.5 * dot(z, u);
  //# priori b0  
  L += - 0.5 * pow( (b0 - mu_b0)/s_b0, 2 );
  //# priori b1
  L += (a_b1 - 1) * log(1 + b1) + (b_b1 - 1) * log(1 - b1);
  //# priori b2
  L += - 0.5 * pow( (b2 - mu_b2)/s_b2, 2 );
  //# jacobiano
  //L += log(1 - pow(b1, 2) );
  
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
  
  vec grad(3), y = param["y_T"], h = param["h"], l = param["l"];
  
  // construindo os vetores u, v e z
  int T = h.n_elem;
  vec z(T), u(T), v(T);
  
  v = l % exp(-h); 
  u = l % exp(- h ) % y.subvec(0, T-1);
  z = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp( h );
  
  grad[0] = dot(v, z) - (b0 - mu_b0) * pow(s_b0, -2);
  
  grad[1] = (1 - pow(b1, 2) ) * ( dot(u, z) ) + (a_b1 - 1) * (1 - b1);
  grad[1] +=  - (b_b1 - 1) * (1 + b1); 
  //# jacobiano
  //grad[1] +=  - 2 * b1;
  
  grad[2] = dot(l, z) - (b2 - mu_b2) * pow(s_b2, -2);
  
  return grad;
}
mat G_b(vec b, List param){
  //b = (b0, delta, b2)
  //param = (y, h, l, mu_b0, s_b0, a_b1, b_b1, mu_b2, s_b2)
  //y = (y0. y1, ..., yT)
  
  double b0 = b[0], b1 = tanh(b[1]), b2 = b[2], s_b0 = param["s_b0"], 
                                                            a_b1 = param["a_b1"], b_b1 = param["b_b1"], s_b2 = param["s_b2"];
  
  vec y = param["y_T"], h = param["h"], l = param["l"];
  int T = h.n_elem;
  
  mat G(3, 3, fill::zeros);
  
  //1a col
  G(0, 0) = sum( l % exp( -h ) ) + pow(s_b0, -2);
  G(1, 0) = (1 - pow(b1, 2) ) * sum( l % exp( -h ) % y.subvec(0, T-1) );
  G(2, 0) = sum( l );
  
  //2a col
  G(0, 1) = G(1, 0);
  G(1, 1) = pow(1 - pow(b1, 2), 2) * sum( l % exp( -h ) % pow(y.subvec(0, T-1), 2) ); 
  G(1, 1) +=  ( 1 - pow(b1, 2) ) * (a_b1 + b_b1 - 2);
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
  
  vec y = param["y_T"], h = param["h"], l = param["l"];
  int T = h.n_elem;
  
  mat dG(3, 3, fill::zeros);
  
  //1a col
  dG(1, 0) = -2 * b1 * (1 - pow(b1, 2) ) * sum( l % exp( -h ) % y.subvec(0, T-1) );
  
  //2a col
  dG(0, 1) = dG(1, 0);
  dG(1, 1) = -2 * b1 * (1 - pow(b1, 2) ) * ( 2 * (1 - pow(b1, 2) ) * sum( l % exp( -h ) % pow(y.subvec(0, T-1), 2) ) + a_b1 + b_b1 - 2);
  dG(2, 1) = -2 * b1 * (1 - pow(b1, 2) ) * sum( l % y.subvec(0, T-1) );   
  
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
//#############################################################################
// ################################ Algorithm #################################
//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);
typedef mat    ( * mat_ptr )(vec, List);
// ############################## hmc function
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
double H_hmc(vec theta, vec p, List param, mat inv_M, num_ptr fun){
  // M = M.i()
  double u, k;
  vec y = param["y_T"];
  
  k = 0.5 * dot(p, inv_M * p);
  u = - fun(theta, param);
  
  return u + k;
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
  
  int acc = 0, dim = theta_init.n_elem; //a = floor(0.1*N)
  double h; //it = N
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
  
  //Progress
  //if( (i % a) == 0 ) cout << "Progresso em " << ceil(100*i/N)<<" %"<< endl;
  
  //double time = timer.toc();
  
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
List Rmhmc(int N, vec eps, int min_L, int max_L, vec theta_init, List param, int fixed_p, 
           num_ptr fun1, num_ptr log_Jac, vec_ptr fun2, mat_ptr M, mat_ptr v[], int param_id){
  
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
    if( M(theta_init, param).is_sympd() ){
      p_current  = mvnrnd(mu, M(theta_init, param) );  
    }else{
      cout << "Error in parameter: " << param_id << "!" << endl;
    }
    
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
List rmhmc_one(vec eps, int min_L, int max_L, vec theta_init, List param, int fixed_p, 
               num_ptr fun1, num_ptr log_Jac, vec_ptr fun2, mat_ptr M, mat_ptr v[], int param_id){
  
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
  if( M(theta_init, param).is_sympd() ){
    p_current  = mvnrnd(mu, M(theta_init, param) );  
  }else{
    cout << "Error in parameter: " << param_id << "!" << endl;
  }
  
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
  
  set_seed( param["seed"] );
  
  /*
  //one draw
  mat_ptr v[3] = {&dG_b_b0, &dG_b_b1, &dG_b_b2};
  List out = rmhmc_one(eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost_b, &log_Jac_b, &glogpost_b, &G_b, v);
  */
  
  //multiple draws
  mat_ptr v[3] = {&dG_b_b0, &dG_b_b1, &dG_b_b2};
  List out = Rmhmc(N, eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost_b, &log_Jac_b, &glogpost_b, &G_b, v, 20);
  
  return out;
  
}
