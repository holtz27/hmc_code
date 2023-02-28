// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double log_gammapdf(double x, double a, double b){
  return log(pow(b, a)*pow(x, (a-1))*exp(-b*x)/tgamma(a));
}
vec U(const vec y, const vec theta, const List param){
  
  double x, mu_0 = param["mu_0"], s_0 = param["s_0"];
  vec Y = log_normpdf(y, theta[0], sqrt(theta[1]));
  
  x = sum(Y);
  x += log_normpdf(theta[0], mu_0, s_0);
  x += log_gammapdf(theta[1], param["a"], param["b"]);
  
  return vec(1, fill::value(-x));
}
vec grad_U(const vec y, const vec theta, const List param){
  
  int N = y.n_elem;
  double grad_mu, grad_sigma2;
  double mu_0=param["mu_0"], s_0=param["s_0"], a=param["a"], b=param["b"];
  
  vec M(N, fill::value(theta[0])), V(N), grad(2);
  
  grad[0] = -sum(y-M)/theta[1] + (theta[0]-mu_0)/pow(s_0, 2);
  
  for(int i = 0; i < N; i++) V[i] = pow(y[i] - M[i], 2);
  
  grad[1] = -0.5*sum(V)/pow(theta[1], 2) + b;
  grad[1] += 0.5*N/theta[1] + (1-a)/theta[1];
  
  return grad;
}

// mu
vec U_mu(const vec y, const vec theta, const List param){
  
  double x, mu_0 = param["mu_0"], s_0 = param["s_0"], sigma2 = param["sigma2"];
  vec Y = log_normpdf(y, theta[0], sqrt(sigma2));
  
  x = sum(Y);
  x += log_normpdf(theta[0], mu_0, s_0);
  //x += log_gammapdf(theta[1], param["a"], param["b"]);
  
  return vec(1, fill::value(-x));
}
vec grad_U_mu(const vec y, const vec theta, const List param){
  
  int N = y.n_elem;
  double mu_0=param["mu_0"], s_0=param["s_0"], sigma2 = param["sigma2"];
  
  vec M(N, fill::value(theta[0])), V(N), grad(1);
  
  grad = -sum(y-M)/sigma2 + (theta[0]-mu_0)/pow(s_0, 2);
  
  //for(int i = 0; i < N; i++) V[i] = pow(y[i] - M[i], 2);
  
  //grad[1] = -0.5*sum(V)/pow(theta[1], 2) + b;
  //grad[1] += 0.5*N/theta[1] + (1-a)/theta[1];
  
  return grad;
}

//sigma2
vec U_sigma2(const vec y, const vec theta, const List param){
  
  double x, mu_0 = param["mu_0"], s_0 = param["s_0"];
  vec Y = log_normpdf(y, theta[0], sqrt(theta[1]));
  
  x = sum(Y);
  //x += log_normpdf(theta[0], mu_0, s_0);
  x += log_gammapdf(theta[1], param["a"], param["b"]);
  
  return vec(1, fill::value(-x));
}
vec grad_U_sigma2(const vec y, const vec theta, const List param){
  
  int N = y.n_elem;
  double grad_mu, grad_sigma2;
  double a = param["a"], b = param["b"];
  
  vec V(N), grad(1);
  
  //grad[0] = -sum(y-M)/theta[1] + (theta[0]-mu_0)/pow(s_0, 2);
  
  for(int i = 0; i < N; i++) V[i] = pow(y[i] - theta[0], 2);
  
  grad[1] = -0.5*sum(V)/pow(theta[1], 2) + b;
  grad[1] += 0.5*N/theta[1] + (1 - a)/theta[1];
  
  return grad;
}

//function pointer
typedef vec (*fun_ptr)(vec, vec, List);

// hmc function
mat lf(const double eps, const int L, const vec y, vec theta_current, 
       const vec p_current, mat M, const List param, fun_ptr fun){
  
  int T = theta_current.n_rows, dim = theta_current.n_elem;
  vec theta = theta_current, p = p_current;
  mat P(dim, 2, fill::zeros);
  
  // Integration Leapfrog
  for(int k = 0; k < L; k++){
    //p -= 0.5 * eps * grad_U(y, theta, param);
    p -= 0.5 * eps * fun(y, theta, param);
    theta += eps * M.i() * p;
    //p -= 0.5 * eps * grad_U(y, theta, param);
    p -= 0.5 * eps * fun(y, theta, param);
  }
  
  P.col(0) = theta;
  P.col(1) = p;
  
  return P;
}

double H(const vec p, const mat M, const vec y, vec theta, 
                const List param, fun_ptr fun){
  double u, k;
  
  k = 0.5*dot(p, M.i() * p);
  //u = as_scalar(U(y, theta, param));
  u = as_scalar(fun(y, theta, param));
  
  return u + k;
}

List hmc_in(const int N, const double eps, const int L, const vec y, 
         vec theta_init, const List param, fun_ptr U, fun_ptr grad_U){
  
  int acc = 0, a = floor(0.1*N), dim = theta_init.n_elem;
  double h, it = N;
  mat chain(dim, N + 1, fill::zeros);
  //inicialindo a cadeia
  chain.col(0) = theta_init;
  vec mean(dim, fill::zeros), V(dim, fill::ones), theta_current = theta_init;
  mat M = diagmat(V), Proposal(dim, 2, fill::ones);
  
  for(int i = 1; i <= N; i++){
    //(i) gerando p_current e u
    double u = R::runif(0, 1);
    vec p_current = mvnrnd(mean, M), 
        p_proposal(dim, fill::zeros), theta_proposal(dim, fill::zeros);
    
    //(ii) Leapfrog 
    Proposal = lf(eps, L, y, theta_current, p_current, M, param, grad_U);
    
    theta_proposal = Proposal.col(0);
    p_proposal = Proposal.col(1);
   
    // (iii) Calculando H's
    h = H(p_current, M, y, theta_current, param, U);
    h -= H(p_proposal, M, y, theta_proposal, param, U);
    
    //(iv) Prop aceitação
    if( u < std::min(1.0, exp(h)) ){
      chain.col(i) = theta_proposal;
      theta_current = theta_proposal;
      acc += 1;
    }else{
      chain.col(i) = theta_current;
    }
    //Progress
    //if( (i % a) == 0 ) cout <<"Progresso em "<<ceil(100*i/N)<<" %"<< endl;
  }
  return List::create(Named("chain") = chain, Named("acc") = acc);
}

// [[Rcpp::export]]
List hmc(const int N, const double eps, const int L, const vec y, 
         vec theta_init, const List param){
  
  vec mu_init = vec(1, fill::value(theta_init[0]));
  vec sigma2_init = vec(1, fill::value(theta_init[1]));
  
  List out = hmc_in(N, eps, L, y, mu_init, param, &U_mu, &grad_U_mu);
  
  return out;
  
}
