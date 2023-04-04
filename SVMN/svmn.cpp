#define M_2PI   6.28318530717958647692528676655
#define M_PI 3.141592653589793238462643383280

#define ZTOL (DBL_EPSILON*10.0)
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>
#include <gsl/gsl_cdf.h>
#include <armadillo>

using namespace std;
using namespace arma;
#include "svmn.h"

int g_dim;
double g_mu0;
double g_sigma2mu;
double g_aphi;
double g_bphi;
double g_asigma;
double g_bsigma;
double g_abeta1;
double g_bbeta1;
double g_beta0;
double g_sigma2beta0;
double g_beta2;
double g_sigma2beta2;
/******************************************************************************/
/******************************************************************************/
const gsl_rng_type * T;
    gsl_rng *r;
unsigned long ss0;

mat svmnsim(vec beta,double mu,double phi,double sigma2,double y0){
mat result=zeros<mat>(g_dim+1,2);


result(0,0)=y0;
result(1,1)=mu+sqrt(sigma2/(1-phi*phi))*gsl_ran_ugaussian(r);
result(1,0)=beta[0]+beta[1]*y0+beta[2]*exp(result(1,1))+exp(0.5*result(1,1))*gsl_ran_ugaussian(r);

int j=1;
for (j=1;j<g_dim+1;j++){

result(j,1)=mu+phi*(result(j-1,1)-mu)+sqrt(sigma2)*gsl_ran_ugaussian(r);
result(j,0)=beta[0]+beta[1]*result(j-1,0)+beta[2]*exp(result(j,1))+exp(0.5*result(j,1))*gsl_ran_ugaussian(r);

}

return result;
}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
vec mvrgaussian(int n){
vec xr=zeros<vec>(n,1);
int j;
for (j=0;j<n;j++){
xr[j]=gsl_ran_ugaussian(r);
}
return xr;}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
double SVMNHMC::lposttheta1(vec theta1, vec vH){
double lptheta1=0.0;
double mu=theta1[0];
double phi=tanh(theta1[1]);
double sigma=exp(theta1[2]);

int j=1;
for (j=1;j<g_dim;j++){
lptheta1+=(-0.5/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu))*(vH[j+1]-mu-phi*(vH[j]-mu));
}

lptheta1+=(0.5*log(1-phi*phi)-(0.5*(1-phi*phi)/(sigma*sigma))*(vH[1]-mu)*(vH[1]-mu));
lptheta1+=((g_aphi-1)*log(1+phi)+(g_bphi-1)*log(1-phi)-(0.5/g_sigma2mu)*(mu-g_mu0)*(mu-g_mu0));
lptheta1+=(-2*(g_asigma+1)*log(sigma)-(g_bsigma/(sigma*sigma))-g_dim*log(sigma));
lptheta1+=(log(1-phi*phi)+log(sigma));

return lptheta1;}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
double SVMNHMC::Gtheta1(vec theta1, mat &G, mat &invG, mat& Gparmu, mat &Gparomega, mat &Gpargamma){
//double mu=theta1[0];
double phi=tanh(theta1[1]);
double sigma=exp(theta1[2]);
G=zeros<mat>(3,3);
invG=zeros<mat>(3,3);
Gparmu=zeros<mat>(3,3);
Gparomega=zeros<mat>(3,3);
Gpargamma=zeros<mat>(3,3);

G(0,0)=((1-phi)*(1-phi)*(g_dim-1)+(1-phi*phi))/(sigma*sigma)+(1.0/g_sigma2mu);
G(1,1)=2*phi*phi+(g_dim-1)*(1-phi*phi)+(g_aphi+g_bphi)*(1-phi*phi);
G(1,2)=2*phi;
G(2,1)=2*phi;
G(2,2)=2*g_dim+(4*g_bsigma/(sigma*sigma));
invG=inv_sympd(G);

Gparomega(0,0)=(-2*(1-phi*phi)/(sigma*sigma))*((g_dim-1)*(1-phi)+phi);
Gparomega(1,1)=2*phi*(1-phi*phi)*(2-g_dim+1-g_aphi-g_bphi);
Gparomega(1,2)=2*(1-phi*phi);
Gparomega(2,1)=2*(1-phi*phi);


Gpargamma(0,0)=-(2/(sigma*sigma))*((g_dim-1)*(1-phi)*(1-phi)+1-phi*phi);
Gpargamma(2,2)=-8*g_bsigma/(sigma*sigma);
return 1;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::gradphitheta1(vec theta1, vec vH,mat invG,mat Gmu,mat Gomega,mat Ggamma){
double mu=theta1[0];
double phi=tanh(theta1[1]);
double sigma=exp(theta1[2]);

vec gradlptheta1=zeros<vec>(3,1);

gradlptheta1[0]=((1-phi*phi)/(sigma*sigma))*(vH[1]-mu)-(mu-g_mu0)/g_sigma2mu;
gradlptheta1[1]=-phi+(phi*(1-phi*phi)/(sigma*sigma))*(vH[1]-mu)*(vH[1]-mu)+(g_aphi-1)*(1-phi)-(g_bphi-1)*(1+phi)-2*phi;
gradlptheta1[2]=-g_dim+((1-phi*phi)/(sigma*sigma))*(vH[1]-mu)*(vH[1]-mu)-2*(g_asigma+1)+2*(g_bsigma/(sigma*sigma))+1;
int j=1;
for (j=1;j<g_dim;j++){
gradlptheta1[0]+=((1-phi)/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu));
gradlptheta1[1]+=((1-phi*phi)/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu))*(vH[j]-mu);
gradlptheta1[2]+=(1/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu))*(vH[j+1]-mu-phi*(vH[j]-mu));
}

gradlptheta1[0]+=(-0.5*trace(invG*Gmu));
gradlptheta1[1]+=(-0.5*trace(invG*Gomega));
gradlptheta1[2]+=(-0.5*trace(invG*Ggamma));
return -gradlptheta1;

}
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::gradHmomtheta1(mat invG,vec p){
return invG  *p;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

vec SVMNHMC::nuHtheta1(mat Gmu,mat Gomega,mat Ggamma,mat invG, vec p){
vec iGp = invG * p;
vec aux1 = iGp.t() * Gmu * iGp;
vec aux2 = iGp.t() * Gomega * iGp;
vec aux3 = iGp.t() * Ggamma * iGp;
vec nuH = zeros<vec>(3,1);
nuH[0] = aux1[0];
nuH[1] = aux2[0];
nuH[2] = aux3[0];
return nuH;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::updatetheta1(vec theta1prev, vec vH, int nfix, int nlpfrog, 
                          double epsilon){
  
vec theta1p = zeros<vec>(3,1);
vec theta1 = zeros<vec>(3,1);
int jfix = 1;
int jgl = 1;
mat invG = zeros<mat>(3,3);
mat Gparmu = zeros<mat>(3,3);
mat Gparomega = zeros<mat>(3,3);
mat Gpargamma = zeros<mat>(3,3);
mat G = zeros<mat>(3,3);
/****************************************************************/
Gtheta1(theta1prev, G, invG, Gparmu, Gparomega, Gpargamma);
vec pcur = chol(G).t() * mvrgaussian(3);
/****************************************************************/
vec aux = pcur.t() * invG * pcur;
vec pa = zeros<vec>(3,1);
vec pb = zeros<vec>(3,1);
vec pc = pcur;
vec gdpth1;
vec gdmom1;
vec gdmom2;
/****************************************************************/
double uMH = 0.0;
double aMH = 0.0;
double phip = 0.0;
double sigmap = 0.0;

double HMtheta1 = -lposttheta1(theta1prev, vH) + 0.5 * log( det(G) ) + 0.5 * aux[0];
double HMthetaf = 0.0;

vec gradHMCth1 = zeros<mat>(3,1);
vec theta1c = theta1prev;
vec theta1a;
vec theta1b;
theta1a = theta1c;
pa = pc;
/**********************************************************/
/**********************************************************/
for( jgl=1;jgl < nlpfrog + 1; jgl++ ){

  gdpth1 = gradphitheta1(theta1a,vH,invG,Gparmu,Gparomega,Gpargamma);
  gradHMCth1 = gdpth1 - 0.5 * nuHtheta1(Gparmu,Gparomega,Gpargamma,invG,pa);
  
/****************************************/
/***** pn+1/2 *****************************/

  for( jfix = 1; jfix < nfix + 1; jfix++){
    pb = pa - 0.5 * epsilon * gradHMCth1;
    //pa=pb;
    gradHMCth1 = gdpth1 - 0.5 * nuHtheta1(Gparmu, Gparomega, Gpargamma, invG, pb);
  }
  
/**********************************************************/
  Gtheta1(theta1a, G, invG, Gparmu, Gparomega, Gpargamma);
  gdmom1 = gradHmomtheta1(invG, pb);
  gdmom2 = gradHmomtheta1(invG, pb);
/****************************************************/
/****thetan+1****************************************/

  for( jfix = 1; jfix < nfix + 1; jfix++){
    
    theta1b = theta1a + 0.5 * epsilon * gdmom1 + 0.5 * epsilon * gdmom2;
    //theta1a=theta1b;
    Gtheta1(theta1b, G, invG, Gparmu, Gparomega, Gpargamma);
    gdmom2 = gradHmomtheta1(invG, pb);
    
  }
/*************************************************************/
  gdpth1 = gradphitheta1(theta1b, vH, invG, Gparmu, Gparomega, Gpargamma);
  gradHMCth1 = gdpth1 - 0.5 * nuHtheta1(Gparmu, Gparomega, Gpargamma, invG, pb);
  pb = pb - 0.5 * epsilon * gradHMCth1;
  pa = pb;
  theta1a = theta1b;
}


aux=pb.t()*invG*pb;

HMthetaf=-lposttheta1(theta1b,vH)+0.5*log(det(G))+0.5*aux[0];
uMH=gsl_rng_uniform(r);

aMH=min(1.0,exp(-HMthetaf+HMtheta1));
cout<<aMH<<endl;
if (uMH<aMH){theta1=theta1b;}
else {theta1=theta1prev;}


/**********************************************************/
/**********************************************************/
return theta1;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
double SVMNHMC::lposttheta2(vec theta2, vec vH){
double lptheta2=0.0;
double beta0=theta2[0];
double beta1=tanh(theta2[1]);
double beta2=theta2[2];
lptheta2+=((-0.5/g_sigma2beta0)*(beta0-g_beta0)*(beta0-g_beta0));
lptheta2+=((g_abeta1-1)*log(1+beta1)+(g_bbeta1-1)*log(1-beta1));
lptheta2+=((-0.5/g_sigma2beta2)*(beta2-g_beta2)*(beta2-g_beta2));
int j=1;
for (j=1;j<g_dim+1;j++){
lptheta2+=((-0.5*exp(-vH[j]))*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j]))*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j])));
}
lptheta2+=log(1-beta1*beta1);
return lptheta2;}
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
double SVMNHMC::Gtheta2(vec theta2,vec vH,mat &G, mat &invG, mat& Gparbeta0,mat &Gpardelta,mat &Gparbeta2){
//double beta0=theta2[0];
double beta1=tanh(theta2[1]);
//double beta2=theta2[2];
G=zeros<mat>(3,3);
invG=zeros<mat>(3,3);
Gparbeta0=zeros<mat>(3,3);
Gpardelta=zeros<mat>(3,3);
Gparbeta2=zeros<mat>(3,3);

G(0,0)=sum(exp(-vH(span(1,g_dim),0)))+(1/g_sigma2beta0);

G(2,2)=sum(exp(vH(span(1,g_dim),0)))+(1/g_sigma2beta2);
G(1,2)=(1-beta1*beta1)*sum(vY(span(0,g_dim-1),0));
G(2,1)=G(1,2);
int j=1;
for(j=1;j<g_dim+1;j++){
G(1,1)+=(1-beta1*beta1)*(1-beta1*beta1)*exp(-vH[j])*vY[j-1]*vY[j-1];
G(0,1)+=(1-beta1*beta1)*exp(-vH[j])*vY[j-1];
}

Gpardelta(1,1)=-4*beta1*G(1,1)-2*beta1*(g_abeta1+g_bbeta1)*(1-beta1*beta1);
G(1,1)+=(g_abeta1+g_bbeta1)*(1-beta1*beta1);
G(1,0)=G(0,1);
G(0,2)=g_dim;
G(2,0)=g_dim;

invG=inv_sympd(G);


Gpardelta(0,1)=-2*beta1*G(0,1);
Gpardelta(1,0)=Gpardelta(0,1);
Gpardelta(1,2)=-2*beta1*G(1,2);
Gpardelta(2,1)=Gpardelta(1,2);
return 1;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::gradphitheta2(vec theta2, vec vH,mat invG,mat Gbeta0,mat Gdelta,mat Gbeta2){
double beta0=theta2[0];
double beta1=tanh(theta2[1]);
double beta2=theta2[2];

vec gradlptheta2=zeros<vec>(3,1);

gradlptheta2[0]=-(beta0-g_beta0)/g_sigma2beta0;
gradlptheta2[1]=(g_abeta1-1)*(1-beta1)-(g_bbeta1-1)*(1+beta1)-2*beta1;
gradlptheta2[2]=-(beta2-g_beta2)/g_sigma2beta2;
int j=1;
for (j=1;j<g_dim+1;j++){
gradlptheta2[0]+=(exp(-vH[j])*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j])));
gradlptheta2[1]+=((1-beta1*beta1)*exp(-vH[j])*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j])))*vY[j-1];
gradlptheta2[2]+=(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j]));
}

gradlptheta2[0]+=(-0.5*trace(invG*Gbeta0));
gradlptheta2[1]+=(-0.5*trace(invG*Gdelta));
gradlptheta2[2]+=(-0.5*trace(invG*Gbeta2));
return -gradlptheta2;

}
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::gradHmomtheta2(mat invG,vec p){
return invG*p;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::nuHtheta2(mat Gbeta0,mat Gdelta,mat Gbeta2,mat invG, vec p){
vec iGp=invG*p;
vec aux1=iGp.t()*Gbeta0*iGp;
vec aux2=iGp.t()*Gdelta*iGp;
vec aux3=iGp.t()*Gbeta2*iGp;
vec nuH=zeros<vec>(3,1);
nuH[0]=aux1[0];
nuH[1]=aux2[0];
nuH[2]=aux3[0];
return nuH;}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::updatetheta2(vec theta2prev,vec vH,int nfix,int nlpfrog,double epsilon){
vec theta2p=zeros<vec>(3,1);
vec theta2=zeros<vec>(3,1);
int jfix=1;
int jgl=1;
mat invG=zeros<mat>(3,3);
mat Gparbeta0=zeros<mat>(3,3);
mat Gpardelta=zeros<mat>(3,3);
mat Gparbeta2=zeros<mat>(3,3);
mat G=zeros<mat>(3,3);
/****************************************************************/
Gtheta2(theta2prev,vH,G,invG,Gparbeta0,Gpardelta,Gparbeta2);
vec pcur=chol(G).t()*mvrgaussian(3);
/****************************************************************/
vec aux=pcur.t()*invG*pcur;
vec pa=zeros<vec>(3,1);
vec pb=zeros<vec>(3,1);
vec pc=pcur;
vec gdpth2;
vec gdmom1;
vec gdmom2;
/****************************************************************/
double uMH=0.0;
double aMH=0.0;



double HMtheta1=-lposttheta2(theta2prev,vH)+0.5*log(det(G))+0.5*aux[0];
double HMthetaf=0.0;

vec gradHMCth2=zeros<mat>(3,1);
vec theta2c=theta2prev;
vec theta2a;
vec theta2b;
theta2a=theta2c;
pa=pc;
/**********************************************************/
/**********************************************************/
for (jgl=1;jgl<nlpfrog+1;jgl++){

gdpth2=gradphitheta2(theta2a,vH,invG,Gparbeta0,Gpardelta,Gparbeta2);
gradHMCth2=gdpth2-0.5*nuHtheta2(Gparbeta0,Gpardelta,Gparbeta2,invG,pa);
/****************************************/
/*****pn+1/2*****************************/
for(jfix=1;jfix<nfix+1;jfix++){
pb=pa-0.5*epsilon*gradHMCth2;
//pa=pb;
gradHMCth2=gdpth2-0.5*nuHtheta2(Gparbeta0,Gpardelta,Gparbeta2,invG,pb);
}
/**********************************************************/
Gtheta2(theta2a,vH,G,invG,Gparbeta0,Gpardelta,Gparbeta2);
gdmom1=gradHmomtheta2(invG,pb);
gdmom2=gradHmomtheta2(invG,pb);
/****************************************************/
/****thetan+1****************************************/
for(jfix=1;jfix<nfix+1;jfix++){
theta2b=theta2a+0.5*epsilon*gdmom1+0.5*epsilon*gdmom2;
//theta1a=theta1b;
Gtheta2(theta2b,vH,G,invG,Gparbeta0,Gpardelta,Gparbeta2);
gdmom2=gradHmomtheta2(invG,pb);
}
/*************************************************************/
gdpth2=gradphitheta2(theta2b,vH,invG,Gparbeta0,Gpardelta,Gparbeta2);
gradHMCth2=gdpth2-0.5*nuHtheta2(Gparbeta0,Gpardelta,Gparbeta2,invG,pb);
pb=pb-0.5*epsilon*gradHMCth2;
pa=pb;
theta2a=theta2b;
}


aux=pb.t()*invG*pb;

HMthetaf=-lposttheta2(theta2b,vH)+0.5*log(det(G))+0.5*aux[0];
uMH=gsl_rng_uniform(r);

aMH=min(1.0,exp(-HMthetaf+HMtheta1));
cout<<aMH<<endl;

if (uMH<aMH){theta2=theta2b;}
else {theta2=theta2prev;}

return theta2;}

/**********************************************************/
/**********************************************************/
/**********************************************************/
/**********************************************************/
double  SVMNHMC::lpostvH(vec vH,vec theta1,vec theta2){
double lpost=0.0;
double mu=theta1[0];
double phi=tanh(theta1[1]);
double sigma=exp(theta1[2]);
double beta0=theta2[0];
double beta1=tanh(theta2[1]);
double beta2=theta2[2];
lpost+=(-0.5*sum(vH(span(1,g_dim),0)));
lpost+=(-0.5*(1-phi*phi)/(sigma*sigma))*(vH[1]-mu)*(vH[1]-mu);
int j=1;
for (j=1;j<g_dim;j++){
lpost+=(-0.5*exp(-vH[j])*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j]))*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j])));
lpost+=((-0.5/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu))*(vH[j+1]-mu-phi*(vH[j]-mu)));
}
lpost+=(-0.5*exp(-vH[g_dim])*(vY[g_dim]-beta0-beta1*vY[g_dim-1]-beta2*exp(vH[g_dim]))*(vY[g_dim]-beta0-beta1*vY[g_dim-1]-beta2*exp(vH[g_dim])));
return lpost;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::gradlpvH(vec vH,vec theta1,vec theta2){
vec gradvH=zeros<vec>(g_dim+1,1);
double mu=theta1[0];
double phi=tanh(theta1[1]);
double sigma=exp(theta1[2]);
double beta0=theta2[0];
double beta1=tanh(theta2[1]);
double beta2=theta2[2];
double auxy=0.0;

auxy=(vY[1]-beta0-beta1*vY[0]-beta2*exp(vH[1]));
gradvH[1]=-0.5+0.5*exp(-vH[1])*auxy*auxy+beta2*auxy-((1-phi*phi)/(sigma*sigma))*(vH[1]-mu)+(phi/(sigma*sigma))*(vH[2]-mu-phi*(vH[1]-mu));

int j=2;
for (j=2;j<g_dim;j++){
auxy=vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j]);
gradvH[j]=-0.5+0.5*exp(-vH[j])*(auxy*auxy)+beta2*auxy+(phi/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu))-(1/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu));
}
auxy=(vY[g_dim]-beta0-beta1*vY[g_dim-1]-beta2*exp(vH[g_dim]));
gradvH[g_dim]=-0.5+0.5*exp(-vH[g_dim])*auxy*auxy+beta2*auxy-(1/(sigma*sigma))*(vH[g_dim]-mu-phi*(vH[g_dim-1]-mu));

return gradvH;}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
vec SVMNHMC::updatevH(vec vHprev,vec theta1,vec theta2,int nlpgrog,double epsilon){
vec vH=zeros<vec>(g_dim+1,1);
vec pcur=zeros<vec>(g_dim+1,1);
pcur(span(1,g_dim),0)=mvrgaussian(g_dim);

double Hm1=0.0;
double Hmf=0.0;
double uMH=0.0;
double aMH=0.0;
vec  auxp=pcur(span(1,g_dim),0).t()*pcur(span(1,g_dim),0);
Hm1=-lpostvH(vHprev,theta1,theta2)+0.5*auxp[0];
int jlp=1;
vec pa=pcur;
vec pb=pcur;
vec vHa=vHprev;
vec vHb=vHprev;
for (jlp=1;jlp<nlpgrog+1;jlp++){
pb=pa+0.5*epsilon*gradlpvH(vHa,theta1,theta2);
//cout<<gradlpvH(vHa,theta1,theta2)<<endl;
//cin>>sk;
vHb=vHa+epsilon*pb;
pb=pb+0.5*epsilon*gradlpvH(vHb,theta1,theta2);
//cout<<+gradlpvH(vHb,theta1,theta2)<<endl;
vHa=vHb;
pa=pb;
}
auxp=pb(span(1,g_dim),0).t()*pb(span(1,g_dim),0);
Hmf=-lpostvH(vHb,theta1,theta2)+0.5*auxp[0];

aMH=min(1.0,exp(-Hmf+Hm1));
cout<<aMH<<endl;
//cin>>sk;
uMH=gsl_rng_uniform(r);
if(uMH<aMH){vH=vHb;}
else {vH=vHprev;}
/*************************************************************************/

return vH;}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int main ()

 {
 srand((unsigned)time( NULL));
 ss0=rand();
  gsl_rng_env_setup();
    T = gsl_rng_ranlxd1;;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, ss0);
/*************************************************************************/
/*************************************************************************/
//g_dim=2000;
vec b0;
//b0<<0.3<<0.3<<-0.25<<endr;
/*mat Dat=svmnsim(b0,0.005,0.98, 0.04, 0.05);*/

//Dat.save("Dat.txt",raw_ascii);
mat Dat;
Dat.load("ibovespa.txt",raw_ascii);
g_dim=Dat.n_rows-1;
g_mu0=0.0;
g_sigma2mu=10;
g_aphi=20.0;
g_bphi=1.5;
g_asigma=10.0;
g_bsigma=0.1;
g_beta0=0.0;
g_beta2=0.0;
g_sigma2beta0=10;
g_sigma2beta2=10;
g_abeta1=5.0;
g_bbeta1=1.5;
/*************************************************************************/
SVMNHMC svmnhmc;
svmnhmc.vY=Dat.col(0);
vec vH=0.0*ones<vec>(g_dim+1,1);;

//vH.load("vHmini.txt",raw_ascii);//=0.05*ones<vec>(g_dim+1,1);//0.7*Dat.col(1);
int burnin=10000;
int fin=20000;
vec theta1=zeros<vec>(3,1);
theta1[0]=0.005;
theta1[1]=0.5*(log(1+0.98)-log(1-0.98));
theta1[2]=log(sqrt(0.017));
vH[1]=0.005+sqrt(0.03)/(1-0.95*0.95)*gsl_ran_ugaussian(r);
int kt=1;
for(kt=2;kt<g_dim+1;kt++){
vH[kt]=0.005+0.95*(vH[kt-1]-0.005)+sqrt(0.03)*gsl_ran_ugaussian(r);
}


int j=1;

vec theta2=zeros<vec>(3,1);
theta2[0]=0.3;
theta2[1]=0.5*(log(1+0.03)-log(1-0.03));
theta2[2]=-0.025;

ofstream theta10file;
ofstream theta20file;
theta10file.open("theta10.txt");
theta20file.open("theta20.txt");

for (j=1;j<burnin+1;j++){
cout<<j<<endl;
theta1=svmnhmc.updatetheta1(theta1,vH,5,20,0.5);
theta2=svmnhmc.updatetheta2(theta2,vH,5,20,0.1);
vH=svmnhmc.updatevH(vH,theta1,theta2,50,0.02);
//vHmean+=vH;

theta10file<<theta1.t();
theta20file<<theta2.t();
}
theta10file.close();
theta20file.close();


ofstream vHfile;
ofstream evHfile;
ofstream e05vHfile;
vHfile.open("outvH.txt");
evHfile.open("outevH.txt");
e05vHfile.open("oute05vH.txt");

ofstream theta11file;
ofstream theta21file;
theta11file.open("theta11.txt");
theta21file.open("theta21.txt");
int rr=0;
vec vHmean=zeros<vec>(g_dim+1,1);
for (j=1+burnin;j<burnin+fin+1;j++){
cout<<j<<endl;
theta1=svmnhmc.updatetheta1(theta1,vH,5,20,0.5);
theta2=svmnhmc.updatetheta2(theta2,vH,5,20,.1);
vH=svmnhmc.updatevH(vH,theta1,theta2,50,0.02);
rr=(j%10);
    if (rr==1){
vHfile<<vH.t();
evHfile<<exp(vH).t();
e05vHfile<<exp(0.5*vH).t();
              }

vHmean+=vH;

theta11file<<theta1.t();
theta21file<<theta2.t();
}
theta11file.close();
theta21file.close();
vHfile.close();
evHfile.close();
e05vHfile.close();

vHmean=(1./fin)*vHmean;
vHmean.save("vHmean.txt",raw_ascii);

    return 1;}




/******************************************************************************/
/*****************************************************************************/
