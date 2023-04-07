#ifndef SVMN_RCPP_BRUNO_H_INCLUDED
#define SVMN_RCPP_BRUNO_H_INCLUDED

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat svmnsim(vec beta, double mu, double phi, double sigma2, double y0);
vec mvrgaussian(int n);

// theta
// Model functions
double logpost_theta(vec theta, vec h, int T);
vec glogpost_theta(vec theta, vec h, mat inv_G, mat dG_mu, mat dG_omega, mat dG_gamma, int T);
// Metrics
double G_theta(vec theta, mat &G, mat &inv_G, mat &dG_mu, mat &dG_omega, mat &dG_gamma, int T);
// Energy functions
vec nuH_theta(mat dG_mu, mat dG_omega, mat dG_gamma, mat inv_G, vec p);
vec gradHmom_theta(mat inv_G, vec p);
// hmc update
vec rmhmc_theta(vec theta_cur, vec h, int fixp, int L, double eps, int T, int &acc);


// beta
// Model functions
double logpost_b(vec b, vec h, int T, vec y_T);
vec glogpost_b(vec b, vec h, mat inv_G, mat G_b0, mat dG_delta, mat dG_b2, int T, vec y_T);
// Metrics
double G_b(vec b, vec h, mat &G, mat &inv_G, mat& dG_b0, mat &dG_delta, mat &dG_b2, int T, vec y_T);
// Energy functions
vec nuH_b(mat dG_b0, mat dG_delta, mat dG_b2, mat inv_G, vec p);
vec gradHmom_b(mat inv_G, vec p);
// hmc update
vec rmhmc_b(vec b_cur, vec h, int fixp, int L, double eps, int T, vec y_T , int &acc);

// h
// Model functions
double logpost_h(vec h, vec theta, vec b, int T, vec y_T);
vec glogpost_h(vec h, vec theta, vec b, int T, vec y_T);
// hmc update
vec hmc_h(vec h_cur, vec theta, vec b, int L, double eps, int T, vec y_T, int &acc);
#endif // SVMN_RCPP_BRUNO_H_INCLUDED
