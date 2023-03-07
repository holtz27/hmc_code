// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//#############################################################################
//########################## l | v ~ Gama(v/2, v/2)

//########### l_t | y, b, h, v ~ Gama( a_t, b_t ), 1 < t < T
//# a_t = (v + 1)/2
//# b_t =  ( exp( -h_t ) * (y_t - b0 - b1 * y_{t-1} - b2 * exp{h_t})**2 + v )/2
//########################## l
vec l_gibbs(List param){
  //param = (y, h, b, v)
  double v = param["v"];
  vec y = param["y"], h = param["h"], b = param["b"];
  int T = h.n_elem;
  
  vec l(T), u(T);
  
  u = exp(-h) % pow(y.subvec(1, T) - b[0] - b[1] * y.subvec(0, T-1) - b[2] * exp(h), 2);
  
  for(int i = 0; i < T; i++){
    l(i) = R::rgamma(0.5 * (v + 1), pow(0.5 * (u[i] + v), -1) );  
  }
  
  return l;
}
//########################## v.1
//######## Transformação: T(v) = e                         ############
//######## e = log(v)
double logpost_v_1(vec e, List param){
  //param = (l, a_v, b_v)
  double v = exp( as_scalar(e) ), L, a_v = param["a_v"], b_v = param["b_v"];
  vec l = param["l"];
  int T = l.n_elem;
  
  L = 0.5 * T * v * log( 0.5 * v  ) - T * log( tgamma( 0.5 * v )  );
  L += 0.5 * v * sum( log(l) - l );
  //# priori
  L += a_v * log( v ) - b_v * v;
  
  return L;
}
vec glogpost_v_1(vec e, List param){
  //param = (l, a_v, b_v)
  double v = exp( as_scalar(e) ), a_v = param["a_v"], b_v = param["b_v"];
  vec l = param["l"], grad(1);
  int T = l.n_elem;
  
  grad(0) = 0.5 * T * v * log(0.5 * v) + 0.5 * T * v - 0.5 * T * v * R::digamma(0.5 * v);
  grad(0) += 0.5 * v * sum( log(l) - l ); 
  grad(0) += a_v - b_v * v;
  
  return grad;
}
mat G_v_1(vec e, List param){
  //#param = (l, a_v, b_v)
  double v = exp( as_scalar(e) ), b_v = param["b_v"];
  vec l = param["l"];
  int T = l.n_elem;
  mat G(1, 1);
  
  G(0, 0) = 0.25 * T * v * v * R::psigamma(0.5 * v, 1);
  G(0, 0) +=  v * b_v - 0.5 * T * v;
  
  return G;
}
mat dG_v_1(vec e, List param){
  //#param = (l, a_v, b_v)
  double v = exp( as_scalar(e) ), b_v = param["b_v"];
  vec l = param["l"];
  int T = l.n_elem;
  mat dG(1, 1);
  
  dG(0, 0) = 0.5 * T * pow(v, 2) * R::psigamma(0.5 * v, 1);
  dG(0, 0) += 0.125 * T * pow(v, 3) * R::psigamma(0.5 * v, 2);
  dG(0, 0) += b_v * v - 0.5 * T * v;
  
  return dG;
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
  
  mat_ptr v[1] = { &dG_v_1 };
  List out = Rmhmc(N, eps, min_L, max_L, theta_init, param, fixed_p,  
                   &logpost_v_1, &glogpost_v_1, &G_v_1, v);
  
  return out;
  
}
##########################################################################################
path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/v_1_rmhmc.cpp'
Rcpp::sourceCpp(path)

N = 5e3
z = rmhmc(N, eps = 0.1, min_L = 30, max_L = 50, 
            fixed_p = 5, 
            theta_init = log( 2 ), 
            param = list(l = l, a_v = 12, b_v = 8) )

z$acc/N
z$time

chain = unlist(z$chain)
sum( is.na(chain) )
chain = exp(chain)

############################### Convergence analysis
# Trace plots
burn = 100
burned = as.matrix(chain[1, -c(1:burn)])

par(mfrow = c(1, 2))
plot(burned, type = 'l')
plot(acf(burned, lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

# Jumps
lags = 5
jumps = seq(1, N - burn, by = 1 + lags)
burned_lag = burned[jumps]

par(mfrow = c(1, 2))
plot(burned_lag, type = 'l')
plot(acf(burned_lag, lag.max = 100, plot = FALSE)[1:100])
par(mfrow = c(1, 1))

par(mfrow = c(1, 3))
plot(burned_lag, type = 'l')
plot(acf(burned_lag, lag.max = 100, plot = FALSE)[1:100])
hist(burned_lag, breaks = 30, xlab = '', main = '')
par(mfrow = c(1, 1))

data = matrix(c( mean(burned_lag), 
                 quantile(burned_lag, probs = c(0.05, 0.975)) )
              ,nrow = 1, byrow = TRUE)
row.names(data) = c('v_hat')
colnames(data) = c('mean', '5%', '97,5%')
data

