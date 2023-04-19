#ifndef SVM_SMN_TS_H_INCLUDED
#define SVM_SMN_TS_H_INCLUDED

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat svmnsim(vec beta, double mu, double phi, double sigma2, double y0);
vec mvrgaussian(int n);
vec gH_mom(mat inv_G, vec p);

// theta
// Model functions
double logpost_theta(vec theta, vec h, int T);
vec glogpost_theta(vec theta, vec h, mat inv_G, mat dG_mu, mat dG_omega, mat dG_gamma, int T);
// Metrics
int G_theta(vec theta, mat &G, mat &inv_G, mat &dG_mu, mat &dG_omega, mat &dG_gamma, int T);
// Energy functions
vec nuH_theta(mat dG_mu, mat dG_omega, mat dG_gamma, mat inv_G, vec p);
//vec gradHmom_theta(mat inv_G, vec p);
// hmc updates
vec hmc_theta(vec theta_cur, vec h, int L, double eps, int T, int &acc);
vec rmhmc_theta(vec theta_cur, vec h, int fixp, int L, vec eps, int T, int &acc);

// beta
// Model functions
double logpost_b(vec b, vec h, vec l, int T, vec y_T);
vec glogpost_b(vec b, vec h, vec l, mat inv_G, mat G_b0, mat dG_delta, mat dG_b2, int T, vec y_T);
// Metrics
int G_b(vec b, vec h, vec l, mat &G, mat &inv_G, mat& dG_b0, mat &dG_delta, mat &dG_b2, int T, vec y_T);
// Energy functions
vec nuH_b(mat dG_b0, mat dG_delta, mat dG_b2, mat inv_G, vec p);
//vec gradHmom_b(mat inv_G, vec p);
// hmc update
vec rmhmc_b(vec b_cur, vec h, vec l, int fixp, int L, double eps, int T, vec y_T , int &acc);

// v
// Model functions
double logpost_v(double e, vec l, int T);
double glogpost_v(double e, vec l, int T);
// Metrics
int G_v(vec e, double &G, double &inv_G, double &dG_v, int T);
// Energy functions
double nuH_v(double dG_v, double inv_G, double p);
double gradHmom_v(double inv_G, double p);
// hmc update
double rmhmc_v(double e_cur, vec l, int fixp, int L, double eps, int T, int &acc);

// l
vec l_gibbs(double e, vec y_T, vec h, vec b, int T);

// h
// Model functions
double logpost_h(vec h, vec theta, vec l, vec b, int T, vec y_T);
vec glogpost_v(vec e, vec l, mat inv_G, mat dG_v, int T);
// hmc update
vec hmc_h(vec h, vec theta, vec b, vec l, int L, double eps, int T, vec y_T, int &acc);

#endif // SVM_SMN_TS_H_INCLUDED
