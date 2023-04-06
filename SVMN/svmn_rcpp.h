#ifndef SVMN_RCPP_H_INCLUDED
#define SVMN_RCPP_H_INCLUDED

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

vec mvrgaussian( int n );

double lposttheta1(vec theta1, vec vH, int g_dim);
double Gtheta1(vec theta1, mat &G, mat &invG, mat &Gparmu, mat &Gparomega, mat &Gpargamma, int g_dim);
vec gradphitheta1(vec theta1, vec vH, mat invG, mat Gmu, mat Gomega, mat Ggamma, int g_dim);
vec nuHtheta1(mat Gmu, mat Gomega, mat Ggamma, mat invG, vec p);
vec updatetheta1(vec theta1prev, vec vH, int nfix, int nlpfrog, double epsilon, int g_dim);
vec gradHmomtheta1(mat invG, vec p);

double Gtheta2(vec theta2,vec vH,mat &G, mat &invG, mat &Gparbeta0,mat &Gpardelta,mat &Gparbeta2, int g_dim, vec vY);
double lposttheta2(vec theta2, vec vH, int g_dim, vec vY);
vec gradphitheta2(vec theta2, vec vH,mat invG,mat Gbeta0,mat Gdelta,mat Gbeta2, int g_dim, vec vY);
vec gradHmomtheta2(mat invG,vec p);
vec nuHtheta2(mat Gbeta0,mat Gdelta,mat Gbeta2,mat invG, vec p);
vec updatetheta2(vec theta2prev,vec vH,int nfix,int nlpfrog,double epsilon, int g_dim, vec vY);

double lpostvH(vec vH,vec theta1,vec theta2, int g_dim, vec vY);
vec gradlpvH(vec vH,vec theta1,vec theta2, int g_dim, vec vY);
vec updatevH(vec vHprev,vec theta1,vec theta2,int nlpgrog,double epsilon, int g_dim, vec vY);

#endif // SVMN_RCPP_H_INCLUDED
