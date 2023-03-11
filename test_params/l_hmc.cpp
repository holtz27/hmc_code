// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//#############################################################################
//########################## l
vec l_gibbs(List param){
  //param = (y, h, b, e)
  // e = log(v)
  double e = param["e"], v = exp( e );
  vec y = param["y_T"], h = param["h"], b = param["b"];
  int T = h.n_elem;
  
  vec l(T), u(T);
  
  u = exp(-h) % pow(y.subvec(1, T) - b[0] - b[1] * y.subvec(0, T-1) - b[2] * exp(h), 2);
  
  for(int i = 0; i < T; i++){
    l(i) = R::rgamma(0.5 * (v + 1), pow(0.5 * (u[i] + v), -1) );  
  }
  
  return l;
  
}
//########################## r
double logpost_r(vec r, List param){
  //#h = (h1, ..., hT)
  //#param = (y, h, b, e)
  vec y = param["y_T"], b = param["b"], h = param["h"];
  double L, b0 = b[0], b1 = b[1], b2 = b[2], e = param["e"], v = exp( e );
  int T = h.n_elem;
  
  vec mu = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp( h );
  vec a = 0.5 * mu % mu % exp( -h ) + 0.5 * v;
  vec l = exp( r );
  
  L = 0.5 * (v - 1) * sum( log( l ) ) - sum( a % l );
  
  return L;
}
vec glogpost_r(vec r, List param){
  //#h = (h1, ..., hT)
  //#param = (y, l, theta, b)
  vec y = param["y_T"], b = param["b"], h = param["h"];
  double L, b0 = b[0], b1 = b[1], b2 = b[2], e = param["e"], v = exp( e );
  int T = h.n_elem;
  
  vec mu = y.subvec(1, T) - b0 - b1 * y.subvec(0, T-1) - b2 * exp( h );
  vec a = 0.5 * mu % mu % exp( -h ) + 0.5 * v;
  vec l = exp( r );
  
  vec grad = 0.5 * ( v - 1 ) - a % l;
  
  return grad;
}
double log_Jac_r(vec r, List param){
  return sum( r );  
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
double H_hmc(vec theta, vec p, List param, mat inv_M, num_ptr fun, 
             num_ptr log_Jac){
  // M = M.i()
  double E;
  vec y = param["y_T"];
  
  E = 0.5 * dot(p, inv_M * p);
  E += - fun(theta, param) - log_Jac(theta, param);
  
  return E;
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
         List param, num_ptr fun1, num_ptr log_Jac, vec_ptr fun2){
  
  wall_clock timer;
  timer.tic();
  
  int acc = 0, a = floor( 0.1 * N ), dim = theta_init.n_elem;
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
    h = H_hmc(theta_current, p_current, param, inv_M, fun1, log_Jac);
    h -= H_hmc(Proposal.col(0), Proposal.col(1), param, inv_M, fun1, log_Jac);
    
    //(iv) Prop aceitação
    if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
      chain.col(i) = Proposal.col(0);
      theta_current = Proposal.col(0);
      acc += 1;
    }else{
      chain.col(i) = theta_current;
    }
    //Progress
    if( (i % a) == 0 ) cout << "Progresso em " << ceil(100*i/N)<<" %"<< endl;
  }
  
  double time = timer.toc();
  
  return List::create(Named("chain") = chain, Named("acc") = acc, Named("time") = time);
}

// ############################## set seed function
void set_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}
// [[Rcpp::export]]
List l_hmc(int N,
             double eps_l, int min_L_l, int max_L_l,
             vec init, List param){
  
  int seed = param["seed"];
  
  if( !( seed == 0) ) set_seed( seed );
  
  List out = hmc(N, eps_l, min_L_l, max_L_l, init, param,
                 &logpost_r, &log_Jac_r, &glogpost_r);
  
  return out;
}
