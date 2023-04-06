#include "svmn_rcpp.h"
// parameters priori

double g_mu0 = 0.0;
double g_sigma2mu = 10;
double g_aphi = 20.0;
double g_bphi = 1.5;
double g_asigma = 10.0;
double g_bsigma = 0.1;
double g_beta0 = 0.0;
double g_beta2 = 0.0;
double g_sigma2beta0 = 10;
double g_sigma2beta2 = 10;
double g_abeta1 = 5.0;
double g_bbeta1 = 1.5;
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
vec mvrgaussian(int n){
    vec xr = zeros<vec>(n,1);
    int j;
    for (j=0;j<n;j++){
        //xr[j]=gsl_ran_ugaussian(r);
        xr[j] = randn( );
    }
    return xr;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
double lposttheta1(vec theta1, vec vH, int g_dim){
    double lptheta1=0.0;
    double mu=theta1[0];
    double phi=tanh(theta1[1]);
    double sigma=exp(theta1[2]);

    for ( int j = 1 ; j < g_dim ; j++ ){
        lptheta1+=(-0.5/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu))*(vH[j]-mu-phi*(vH[j-1]-mu));
    }

    lptheta1+=(0.5*log(1-phi*phi)-(0.5*(1-phi*phi)/(sigma*sigma))*(vH[0]-mu)*(vH[0]-mu));
    lptheta1+=((g_aphi-1)*log(1+phi)+(g_bphi-1)*log(1-phi)-(0.5/g_sigma2mu)*(mu-g_mu0)*(mu-g_mu0));
    lptheta1+=(-2*(g_asigma+1)*log(sigma)-(g_bsigma/(sigma*sigma))-g_dim*log(sigma));
    // log Jacobiano
    lptheta1 += ( log(1-phi*phi)+log(sigma));

    return lptheta1;
}
double Gtheta1(vec theta1, mat &G, mat &invG, mat& Gparmu, mat &Gparomega, mat &Gpargamma, int g_dim){
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
    invG = inv_sympd(G);

    Gparomega(0,0)=(-2*(1-phi*phi)/(sigma*sigma))*((g_dim-1)*(1-phi)+phi);
    Gparomega(1,1)=2*phi*(1-phi*phi)*(2-g_dim+1-g_aphi-g_bphi);
    Gparomega(1,2)=2*(1-phi*phi);
    Gparomega(2,1)=2*(1-phi*phi);

    Gpargamma(0,0)=-(2/(sigma*sigma))*((g_dim-1)*(1-phi)*(1-phi)+1-phi*phi);
    Gpargamma(2,2)=-8*g_bsigma/(sigma*sigma);
    return 1;
}
vec gradphitheta1(vec theta1, vec vH, mat invG, mat Gmu, mat Gomega, mat Ggamma, int g_dim){
    double mu=theta1[0];
    double phi=tanh(theta1[1]);
    double sigma=exp(theta1[2]);

    vec gradlptheta1=zeros<vec>(3,1);

    gradlptheta1[0]=((1-phi*phi)/(sigma*sigma))*(vH[0]-mu)-(mu-g_mu0)/g_sigma2mu;
    gradlptheta1[1]=-phi+(phi*(1-phi*phi)/(sigma*sigma))*(vH[0]-mu)*(vH[0]-mu)+(g_aphi-1)*(1-phi)-(g_bphi-1)*(1+phi) - 2*phi;
    gradlptheta1[2]=-g_dim+((1-phi*phi)/(sigma*sigma))*(vH[0]-mu)*(vH[0]-mu)-2*(g_asigma+1)+2*(g_bsigma/(sigma*sigma)) + 1;
    
    for ( int j = 1 ; j < g_dim ; j++ ){
        gradlptheta1[0]+=((1-phi)/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu));
        gradlptheta1[1]+=((1-phi*phi)/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu))*(vH[j-1]-mu);
        gradlptheta1[2]+=(1/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu))*(vH[j]-mu-phi*(vH[j-1]-mu));
    }

    gradlptheta1[0]+=(-0.5*trace(invG*Gmu));
    gradlptheta1[1]+=(-0.5*trace(invG*Gomega));
    gradlptheta1[2]+=(-0.5*trace(invG*Ggamma));
    return -gradlptheta1;

}

vec gradHmomtheta1(mat invG, vec p){
    return invG * p;
}
vec nuHtheta1(mat Gmu, mat Gomega, mat Ggamma, mat invG, vec p){
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
vec updatetheta1(vec theta1prev, vec vH, int nfix, int nlpfrog, double epsilon, int g_dim){
  
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
    Gtheta1(theta1prev, G, invG, Gparmu, Gparomega, Gpargamma, g_dim);
    vec pcur = chol(G).t() * mvrgaussian(3);
    /****************************************************************/
    vec aux = pcur.t() * invG * pcur;
    vec pa = zeros<vec>(3,1);
    vec pb = zeros<vec>(3,1);
    vec pc = pcur;
    vec gdpth1(3);
    vec gdmom1(3);
    vec gdmom2(3);
    /****************************************************************/
    double uMH = 0.0;
    double aMH = 0.0;
    //double phip = 0.0;
    //double sigmap = 0.0;

    //double HMtheta1 = - lposttheta1(theta1prev, vH, g_dim) + 0.5 * log( det(G) ) + 0.5 * aux[0];
    double HMtheta1 = - lposttheta1(theta1prev, vH, g_dim) + 0.5 * log( det(G) ) + 0.5 * as_scalar( aux );
    double HMthetaf = 0.0;

    vec gradHMCth1 = zeros<mat>(3,1);
    vec theta1c = theta1prev;
    vec theta1a;
    vec theta1b;
    theta1a = theta1c;
    pa = pc;
    /**********************************************************/
    /**********************************************************/
    for( jgl = 1;jgl < nlpfrog + 1; jgl++ ){

    //gdpth1 = gradphitheta1(theta1a,vH,invG,Gparmu,Gparomega,Gpargamma, g_dim);
    gradHMCth1 = gdpth1 - 0.5 * nuHtheta1(Gparmu,Gparomega,Gpargamma,invG,pa);
    
    /****************************************/
    /***** pn+1/2 *****************************/

    for( jfix = 1; jfix < nfix + 1; jfix++){
        pb = pa - 0.5 * epsilon * gradHMCth1;
        //pa=pb;
        gradHMCth1 = gdpth1 - 0.5 * nuHtheta1(Gparmu, Gparomega, Gpargamma, invG, pb);
    }
    
    /**********************************************************/
    Gtheta1(theta1a, G, invG, Gparmu, Gparomega, Gpargamma, g_dim);
    gdmom1 = gradHmomtheta1(invG, pb);
    gdmom2 = gradHmomtheta1(invG, pb);
    /****************************************************/
    /****thetan+1****************************************/

    for( jfix = 1; jfix < nfix + 1; jfix++){
        
        theta1b = theta1a + 0.5 * epsilon * gdmom1 + 0.5 * epsilon * gdmom2;
        //theta1a=theta1b;
        Gtheta1(theta1b, G, invG, Gparmu, Gparomega, Gpargamma, g_dim);
        gdmom2 = gradHmomtheta1(invG, pb);
        
    }
    /*************************************************************/
    //gdpth1 = gradphitheta1(theta1b, vH, invG, Gparmu, Gparomega, Gpargamma, g_dim);
    gradHMCth1 = gdpth1 - 0.5 * nuHtheta1(Gparmu, Gparomega, Gpargamma, invG, pb);
    pb = pb - 0.5 * epsilon * gradHMCth1;
    pa = pb;
    theta1a = theta1b;
    }

    aux = pb.t() * invG * pb;

    //HMthetaf = -lposttheta1(theta1b, vH, g_dim) + 0.5 * log( det(G) ) + 0.5 * aux[0];
    HMthetaf = -lposttheta1(theta1b, vH, g_dim) + 0.5 * log( det(G) ) + 0.5 * as_scalar( aux );
    //uMH=gsl_rng_uniform(r);
    uMH = randu( distr_param( 0, 1 ) );

    aMH = std::min( 1.0, exp( -HMthetaf + HMtheta1 ) );
    //cout<<aMH<<endl;
    if ( uMH < aMH ){ 
        theta1 = theta1b;
        }
    else{ 
        theta1 = theta1prev;
        }

    /**********************************************************/
    /**********************************************************/
    return theta1;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
double lposttheta2(vec theta2, vec vH, int g_dim, vec vY){
    double lptheta2=0.0;
    double beta0=theta2[0];
    double beta1=tanh(theta2[1]);
    double beta2=theta2[2];
    lptheta2+=((-0.5/g_sigma2beta0)*(beta0-g_beta0)*(beta0-g_beta0));
    lptheta2+=((g_abeta1-1)*log(1+beta1)+(g_bbeta1-1)*log(1-beta1));
    lptheta2+=((-0.5/g_sigma2beta2)*(beta2-g_beta2)*(beta2-g_beta2));
    
    for ( int j = 1 ; j < g_dim + 1 ; j++ ){
        lptheta2+=((-0.5*exp(-vH[j-1]))*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j-1]))*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j-1])));
    }

    // log Jacobiano
    lptheta2+=log(1-beta1*beta1);
    return lptheta2;
}
double Gtheta2(vec theta2, vec vH, mat &G, mat &invG, mat& Gparbeta0, mat &Gpardelta,mat &Gparbeta2, int g_dim, vec vY){
   //double beta0=theta2[0];
    double beta1=tanh(theta2[1]);
    //double beta2=theta2[2];
    G=zeros<mat>(3,3);
    invG=zeros<mat>(3,3);
    Gparbeta0=zeros<mat>(3,3);
    Gpardelta=zeros<mat>(3,3);
    Gparbeta2=zeros<mat>(3,3);

    G(0,0)=sum(exp(-vH(span(1, g_dim),0)))+(1/g_sigma2beta0);
    G(2,2)=sum(exp(vH(span(1,g_dim),0)))+(1/g_sigma2beta2);
    G(1,2)=(1-beta1*beta1)*sum(vY(span(0,g_dim-1),0));
    G(2,1)=G(1,2);
    
    for( int j = 1 ; j < g_dim + 1 ; j++ ){
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
vec gradphitheta2(vec theta2, vec vH, mat invG, mat Gbeta0, mat Gdelta, mat Gbeta2, int g_dim, vec vY){
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

vec gradHmomtheta2(mat invG, vec p){
    return invG*p;
}
vec nuHtheta2(mat Gbeta0, mat Gdelta, mat Gbeta2, mat invG, vec p){
    vec iGp=invG*p;
    vec aux1=iGp.t()*Gbeta0*iGp;
    vec aux2=iGp.t()*Gdelta*iGp;
    vec aux3=iGp.t()*Gbeta2*iGp;
    vec nuH=zeros<vec>(3,1);
    nuH[0]=aux1[0];
    nuH[1]=aux2[0];
    nuH[2]=aux3[0];
    return nuH;
}
vec updatetheta2(vec theta2prev,vec vH, int nfix, int nlpfrog, double epsilon, int g_dim, vec vY){
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
    Gtheta2(theta2prev, vH, G, invG, Gparbeta0, Gpardelta, Gparbeta2, g_dim, vY);
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

    double HMtheta1=-lposttheta2(theta2prev, vH, g_dim, vY)+0.5*log(det(G))+0.5*aux[0];
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

    gdpth2=gradphitheta2(theta2a, vH, invG, Gparbeta0, Gpardelta, Gparbeta2, g_dim, vY);
    gradHMCth2=gdpth2-0.5*nuHtheta2(Gparbeta0,Gpardelta,Gparbeta2,invG,pa);
    /****************************************/
    /*****pn+1/2*****************************/
    for(jfix=1;jfix<nfix+1;jfix++){
    pb=pa-0.5*epsilon*gradHMCth2;
    //pa=pb;
    gradHMCth2=gdpth2-0.5*nuHtheta2(Gparbeta0,Gpardelta,Gparbeta2,invG,pb);
    }
    /**********************************************************/
    Gtheta2(theta2a,vH,G,invG,Gparbeta0,Gpardelta,Gparbeta2, g_dim, vY);
    gdmom1=gradHmomtheta2(invG,pb);
    gdmom2=gradHmomtheta2(invG,pb);
    /****************************************************/
    /****thetan+1****************************************/
    for(jfix=1;jfix<nfix+1;jfix++){
    theta2b=theta2a+0.5*epsilon*gdmom1+0.5*epsilon*gdmom2;
    //theta1a=theta1b;
    Gtheta2(theta2b,vH,G,invG,Gparbeta0,Gpardelta,Gparbeta2, g_dim, vY);
    gdmom2=gradHmomtheta2(invG,pb);
    }
    /*************************************************************/
    gdpth2=gradphitheta2(theta2b,vH,invG,Gparbeta0,Gpardelta,Gparbeta2, g_dim, vY);
    gradHMCth2=gdpth2-0.5*nuHtheta2(Gparbeta0,Gpardelta,Gparbeta2,invG,pb);
    pb=pb-0.5*epsilon*gradHMCth2;
    pa=pb;
    theta2a=theta2b;
    }

    aux=pb.t()*invG*pb;

    HMthetaf=-lposttheta2(theta2b,vH, g_dim, vY)+0.5*log(det(G))+0.5*aux[0];
    //uMH=gsl_rng_uniform(r);
    uMH = randu( distr_param( 0, 1 ) );

    aMH=std::min(1.0,exp(-HMthetaf+HMtheta1));
    //cout<<aMH<<endl;

    if (uMH<aMH){theta2=theta2b;}
    else {theta2=theta2prev;}

    return theta2;
}
/**********************************************************/
/**********************************************************/
/**********************************************************/
double lpostvH(vec vH, vec theta1, vec theta2, int g_dim, vec vY){
    double lpost=0.0;

    double mu=theta1[0];
    double phi=tanh(theta1[1]);
    double sigma=exp(theta1[2]);
    double beta0=theta2[0];
    double beta1=tanh(theta2[1]);
    double beta2=theta2[2];

    lpost+=(-0.5 * sum( vH( span(0, g_dim-1), 0) ) );
    lpost+=(-0.5*(1-phi*phi)/(sigma*sigma))*(vH[0]-mu)*(vH[0]-mu);
    //int j=1;
    for ( int j = 1 ; j < g_dim ; j++ ){
        lpost+=(-0.5*exp(-vH[j-1])*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j-1]))*(vY[j]-beta0-beta1*vY[j-1]-beta2*exp(vH[j-1])));
        lpost+=((-0.5/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu))*(vH[j]-mu-phi*(vH[j-1]-mu)));
    }
    lpost+=(-0.5*exp(-vH[g_dim-1])*(vY[g_dim]-beta0-beta1*vY[g_dim-1]-beta2*exp(vH[g_dim - 1]))*(vY[g_dim]-beta0-beta1*vY[g_dim-1]-beta2*exp(vH[g_dim-1])));
    return lpost;
}
vec gradlpvH(vec vH, vec theta1, vec theta2, int g_dim, vec vY){
    
    vec gradvH=zeros<vec>(g_dim,1);

    double mu=theta1[0];
    double phi=tanh(theta1[1]);
    double sigma=exp(theta1[2]);
    double beta0=theta2[0];
    double beta1=tanh(theta2[1]);
    double beta2=theta2[2];
    double auxy=0.0;

    auxy = ( vY[1]-beta0-beta1*vY[0]-beta2*exp(vH[0]) );
    gradvH[0]=-0.5+0.5*exp(-vH[0])*auxy*auxy+beta2*auxy-((1-phi*phi)/(sigma*sigma))*(vH[0]-mu)+(phi/(sigma*sigma))*(vH[1]-mu-phi*(vH[0]-mu));

    //int j=2;
    for ( int j = 1 ; j < g_dim - 1 ; j++ ){
        auxy = vY[j+1]-beta0-beta1*vY[j]-beta2*exp(vH[j]);
        gradvH[j] = -0.5 + 0.5*exp(-vH[j])*(auxy*auxy) + beta2 * auxy + (phi/(sigma*sigma))*(vH[j+1]-mu-phi*(vH[j]-mu))-(1/(sigma*sigma))*(vH[j]-mu-phi*(vH[j-1]-mu));
    }

    auxy=(vY[g_dim]-beta0-beta1*vY[g_dim-1]-beta2*exp(vH[g_dim-1]));
    gradvH[g_dim-1]=-0.5+0.5*exp(-vH[g_dim-1])*auxy*auxy  + beta2 * auxy - (1/(sigma*sigma))*(vH[g_dim-1]-mu-phi*(vH[g_dim-2]-mu));

    return gradvH;
}
vec updatevH(vec vHprev,vec theta1,vec theta2,int nlpgrog,double epsilon, int g_dim, vec vY){
    vec vH=zeros<vec>(g_dim+1,1);
    vec pcur=zeros<vec>(g_dim+1,1);
    pcur(span(1,g_dim),0)=mvrgaussian(g_dim);

    double Hm1=0.0;
    double Hmf=0.0;
    double uMH=0.0;
    double aMH=0.0;
    vec  auxp=pcur(span(1,g_dim),0).t()*pcur(span(1,g_dim),0);
    Hm1=-lpostvH(vHprev,theta1,theta2, g_dim, vY)+0.5*auxp[0];
    int jlp=1;
    vec pa=pcur;
    vec pb=pcur;
    vec vHa=vHprev;
    vec vHb=vHprev;
    for (jlp=1;jlp<nlpgrog+1;jlp++){
    pb=pa+0.5*epsilon*gradlpvH(vHa,theta1,theta2, g_dim, vY);
    //cout<<gradlpvH(vHa,theta1,theta2)<<endl;
    //cin>>sk;
    vHb=vHa+epsilon*pb;
    pb=pb+0.5*epsilon*gradlpvH(vHb,theta1,theta2, g_dim, vY);
    //cout<<+gradlpvH(vHb,theta1,theta2)<<endl;
    vHa=vHb;
    pa=pb;
    }
    auxp=pb(span(1,g_dim),0).t()*pb(span(1,g_dim),0);
    Hmf=-lpostvH(vHb,theta1,theta2, g_dim, vY)+0.5*auxp[0];

    aMH=std::min(1.0,exp(-Hmf+Hm1));
    //cout<<aMH<<endl;
    //cin>>sk;
    //uMH=gsl_rng_uniform(r);
    uMH = randu( distr_param( 0, 1 ) );
    if(uMH<aMH){vH=vHb;}
    else {vH=vHprev;}
    /*************************************************************************/

    return vH;
}
