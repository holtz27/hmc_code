// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
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

//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);

// hmc function
mat lf(double eps, int L, vec theta_current, vec p_current, mat inv_M, 
       List param, vec_ptr fun){
  // M = M.i()
  int T = theta_current.n_elem;
  vec theta = theta_current, p = p_current, y = param["y"];
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

double H(vec theta, vec p, List param, mat inv_M, num_ptr fun){
  // M = M.i()
  double u, k;
  vec y = param["y"];
  
  k = 0.5 * dot(p, inv_M * p);
  u = - fun(theta, param);
  
  return u + k;
}

List hmc_in(int N, double eps, int min_L, int max_L, vec theta_init, List param, 
            num_ptr fun1, vec_ptr fun2){
  
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
    h = H(theta_current, p_current, param, inv_M, fun1);
    h -= H(Proposal.col(0), Proposal.col(1), param, inv_M, fun1);
    
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
// [[Rcpp::export]]

List hmc(int N, double eps, int min_L, int max_L, vec theta_init, List param){
  
 List out = hmc_in(N, eps, min_L, max_L, theta_init, param, 
                   &logpost_h, &glogpost_h);
 return out;
  
}
######################################################################################
path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/hmc/h_hmc.cpp'
Rcpp::sourceCpp(path)

N = 3e4
z = hmc(N, eps = 0.015, min_L = 50, max_L = 80,
        theta_init = rep(0, T),
        param = list(y = c(y0, y),
                     l = l, 
                     theta = c(mu, phi, sigma), 
                     b = c(b0, b1, b2) )
        )

z$acc/N
z$time

chain = unlist(z$chain)
sum( is.na(chain) )

####################### Visualization
chain = chain[, 2:(T+1)]
h_hat = apply(chain, MARGIN = 2, FUN = mean)
h_min = apply(chain, MARGIN = 2, FUN = quantile, probs = c(0.05) )
h_max = apply(chain, MARGIN = 2, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:T, h, h_hat, h_min, h_max), ncol = 5)
data = data.frame(data)
names(data) = c('obs', 'vdd','media', 'min','max')
data

library(ggplot2)
g = ggplot(data) 
g = g + geom_line(aes(obs, media))
g = g + geom_line(aes(obs, vdd), color = 'red')
g = g + geom_line(aes(obs, min), linetype = 'dashed')
g = g + geom_line(aes(obs, max), linetype = 'dashed')
g