#ifndef HMC_H
#define HMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);
typedef mat    ( * mat_ptr )(vec, List);

//###########################################################################################################################################################
//###################################################################### HMC ################################################################################
//###########################################################################################################################################################

// Energy function
double H_hmc(vec theta, vec p, List param, mat inv_M, num_ptr fun);

// Leapfrog function
mat lf(double eps, int L, vec theta_current, vec p_current, mat inv_M, List param, vec_ptr fun);

// hmc function
List hmc_one(double eps, int min_L, int max_L, vec theta_init, List param, num_ptr fun1, vec_ptr fun2);

//###########################################################################################################################################################
//#################################################################### RMHMC ################################################################################
//###########################################################################################################################################################

// Energy function 
double H_rmhmc(vec theta, vec p, List param, num_ptr fun, mat_ptr M);

// Grads Energy function
vec dH_p(vec theta, vec p, List param, mat_ptr M);
vec dH_theta(vec theta, vec p, List param, vec_ptr fun, mat_ptr M, mat_ptr v[]);

// Generalized Leapfrog function
mat glf(vec eps, int L, vec theta_current, vec p_current, int fixed_p, List param, vec_ptr fun, mat_ptr M, mat_ptr v[]);

// rmhmc function
List rmhmc_one(vec eps, int min_L, int max_L, vec theta_init, List param, int fixed_p, num_ptr fun1, vec_ptr fun2, mat_ptr M, mat_ptr v[]);


#endif // HMC_H
