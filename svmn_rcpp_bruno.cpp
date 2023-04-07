#include "svmn_rcpp_bruno.h"

// parameters priori
double mu_0 = 0.0, s_0 = 3.2;
double a_phi = 20.0, b_phi = 1.5;
double a_s = 10.0, b_s = 0.1;
double mu_b0 = 0.0, s_b0 = 3.2;
double a_b1 = 10, b_b1 = 10;
double mu_b2 = 5.0, s_b2 = 1.5;

//#############################################################################
//#############################################################################

vec mvrgaussian(int n){
    vec xr = zeros<vec>(n, 1);
    for( int j = 0 ; j < n ; j++ ){
        xr[j] = randn( );
    }
    return xr;
}

//########################## theta = (mu, phi, sigma) #########################
//######## Transformação: T(theta) = theta'                         ###########  
//######## theta' = (mu, arctanh(phi), log(sigma)) = (mu, w, gama) ############

double logpost_theta(vec theta, vec h, int T){
  //theta = (mu, w, gama)
  double L = 0.0; 
  double mu = theta[0];
  double phi = tanh( theta[1] ); 
  double sigma = exp( theta[2] );
  
  L += 0.5 * log(1 - phi * phi ) - T * log( sigma ); 
  L -= 0.5 * (1 - phi * phi ) * ( (h[0] - mu)/sigma ) * ( (h[0] - mu)/sigma ); 
  
  for( int j = 1 ; j < T ; j++ ){
    L += - (0.5 / (sigma * sigma)) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
    
  //# priori mu
  L += - 0.5 * ((mu - mu_0)/s_0) * ((mu - mu_0)/s_0);
  //# priori phi
  L += (a_phi - 1) * log(1 + phi) + (b_phi - 1) * log(1 - phi);
  //# priori sigma
  L += - 2 * (a_s + 1) * log(sigma) - b_s / ( sigma * sigma);
  //# jacobiano de T
  L += log( 1 - phi * phi ) + log( sigma );
  
  return L; 
}
vec glogpost_theta(vec theta, vec h, mat inv_G, mat dG_mu, mat dG_omega, mat dG_gamma, int T){
  //theta = (mu, w, gama)
  double mu = theta[0];
  double phi = tanh( theta[1] ); 
  double sigma = exp( theta[2] );
   
  vec grad = zeros<vec>(3, 1); 
  
  // gradiente mu
  grad[0] += (1 - phi * phi ) * (h[0] - mu) / (sigma * sigma); 
  for( int j = 1 ; j < T ; j++ ){
    grad[0] += (1 - phi) / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
  
  //priori
  grad[0] += - (mu - mu_0) / (s_0 * s_0);
  
  //gradiente w
  grad[1] += - phi + phi * (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma);
  for( int j = 1 ; j < T ; j++ ){
    grad[1] += (1 - phi * phi) / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) ) * (h[j-1] - mu);
  }
  //priori
  grad[1] += (a_phi - 1) * (1 - phi);
  grad[1] -= (b_phi - 1) * (1 + phi);
  //jacobiano
  grad[1] += - 2 * phi;
  
  // gradiente gama
  grad[2] += - T + (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma); 
  for( int j = 1 ; j < T ; j++ ){
    grad[2] += 1.0 / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
  // priori
  grad[2] += - 2 * (a_s + 1) + 2 * b_s / (sigma * sigma);
  // jacobiano
  grad[2] += 1;
  
  // Partial Energy function
  grad[0] += ( -0.5 * trace( inv_G * dG_mu ) );
  grad[1] += ( -0.5 * trace( inv_G * dG_omega ) );
  grad[2] += ( -0.5 * trace( inv_G * dG_gamma ) );
  
  return - grad;
}

double G_theta(vec theta, mat &G, mat &inv_G, mat& dG_mu, mat &dG_omega, mat &dG_gamma, int T){
    
    //theta = (mu, w, gama)
    double phi = tanh( theta[1] );
    double sigma = exp( theta[2] );
    
    G = zeros<mat>(3, 3);
    inv_G = zeros<mat>(3, 3);
    dG_mu = zeros<mat>(3, 3);
    dG_omega = zeros<mat>(3, 3);
    dG_gamma = zeros<mat>(3, 3);

    G(0, 0) = ( (1 - phi ) * ( 1 - phi ) * ( T - 1 ) + ( 1 - phi * phi ) ) / ( sigma*sigma ) + ( 1.0 / (s_0 * s_0) );
    G(1, 1) = 2 * phi * phi + ( T - 1 ) * ( 1 - phi * phi ) + ( a_phi + b_phi ) * ( 1 - phi * phi );
    G(1, 2) = 2 * phi;
    G(2, 1) = 2 * phi;
    G(2, 2) = 2 * T + ( 4 * b_s / ( sigma * sigma ) );
    inv_G = inv_sympd( G );

    dG_omega(0, 0) = ( -2 * ( 1 - phi * phi ) / ( sigma * sigma ) ) * ( ( T - 1 ) * ( 1 - phi ) + phi );
    dG_omega(1, 1) = 2 * phi * ( 1 - phi * phi ) * ( 2 - T + 1 - a_phi - b_phi );
    dG_omega(1, 2) =  2 * ( 1 - phi * phi );
    dG_omega(2, 1) = 2 * ( 1 - phi * phi );

    dG_gamma(0, 0) = - ( 2 / ( sigma * sigma ) ) * ( ( T - 1 ) * ( 1 - phi ) * ( 1 - phi ) + 1 - phi * phi );
    dG_gamma(2, 2) = - 8 * b_s / ( sigma * sigma );
    return 1;
}

vec gradHmom_theta(mat inv_G, vec p){
    return inv_G * p;
}

vec nuH_theta(mat dG_mu, mat dG_omega, mat dG_gamma, mat inv_G, vec p){
    vec iGp = inv_G * p;
    vec aux1 = iGp.t() * dG_mu * iGp;
    vec aux2 = iGp.t() * dG_omega * iGp;
    vec aux3 = iGp.t() * dG_gamma * iGp;
    vec nuH = zeros<vec>(3, 1);
    nuH[0] = aux1[0];
    nuH[1] = aux2[0];
    nuH[2] = aux3[0];
    
    return nuH;
}

vec rmhmc_theta(vec theta_cur, vec h, int fixp, int L, double eps, int T, int &acc){
  
    vec theta1p = zeros<vec>(3, 1);
    vec theta1 = zeros<vec>(3, 1);
    mat inv_G = zeros<mat>(3, 3);
    mat dG_mu = zeros<mat>(3, 3);
    mat dG_omega = zeros<mat>(3, 3);
    mat dG_gamma = zeros<mat>(3, 3);
    mat G = zeros<mat>(3, 3);
    //################################################################
    G_theta(theta_cur, G, inv_G, dG_mu, dG_omega, dG_gamma, T);
    vec pcur = chol( G ).t() * mvrgaussian( 3 );
    //################################################################
    vec aux = pcur.t() * inv_G * pcur;
    vec pa = zeros<vec>(3, 1);
    vec pb = zeros<vec>(3, 1);
    vec pc = pcur;
    vec gdpth1( 3 );
    vec gdmom1( 3 );
    vec gdmom2( 3 );
    //################################################################
    double uMH = 0.0;
    double aMH = 0.0;
    double HMtheta1 = - logpost_theta( theta_cur, h, T ) + 0.5 * log( det(G) ) + 0.5 * aux[0];
    double HMthetaf = 0.0;

    vec gradHMCth1 = zeros<mat>(3, 1);
    vec theta1c = theta_cur;
    vec theta1a;
    vec theta1b;
    theta1a = theta1c;
    pa = pc;
    //################################################################
    //################################################################
    // Generalized Leapfrog function
	
    for( int jgl = 1 ; jgl < L + 1 ; jgl++ ){
	
	gdpth1 = glogpost_theta( theta1a, h, inv_G, dG_mu, dG_omega, dG_gamma, T);
    	gradHMCth1 = gdpth1 - 0.5 * nuH_theta( dG_mu, dG_omega, dG_gamma, inv_G, pa );
    
    	//######################################################
    	//###################### pn + 1/2 ######################
	
    	for( int jfix = 1; jfix < fixp + 1 ; jfix++ ){
        	pb = pa - 0.5 * eps * gradHMCth1;
        	gradHMCth1 = gdpth1 - 0.5 * nuH_theta( dG_mu, dG_omega, dG_gamma, inv_G, pb );
    	}
    
    	//######################################################
	//######################################################
    	G_theta( theta1a, G, inv_G, dG_mu, dG_omega, dG_gamma, T );
    	gdmom1 = gradHmom_theta( inv_G, pb);
    	gdmom2 = gradHmom_theta( inv_G, pb);
    	//######################################################
    	//##################### theta + 1 ######################

    	for( int jfix = 1; jfix < fixp + 1; jfix ++ ){
        
        	theta1b = theta1a + 0.5 * eps * gdmom1 + 0.5 * eps * gdmom2;
        	G_theta( theta1b, G, inv_G, dG_mu, dG_omega, dG_gamma, T );
        	gdmom2 = gradHmom_theta( inv_G, pb );
        
    	}
    	//######################################################
    	gdpth1 = glogpost_theta( theta1b, h, inv_G, dG_mu, dG_omega, dG_gamma, T);
    	gradHMCth1 = gdpth1 - 0.5 * nuH_theta( dG_mu, dG_omega, dG_gamma, inv_G, pb );
    	pb = pb - 0.5 * eps * gradHMCth1;
    	pa = pb;
    	theta1a = theta1b;
    }

    aux = pb.t() * inv_G * pb;

    HMthetaf = -logpost_theta( theta1b, h, T ) + 0.5 * log( det( G ) ) + 0.5 * aux[0];
   
    uMH = randu( distr_param( 0, 1 ) );

    aMH = std::min( 1.0, exp( -HMthetaf + HMtheta1 ) );
    
    if ( uMH < aMH ){ 
        theta1 = theta1b;
        acc ++; 
        }else{ 
        theta1 = theta_cur;
        }

    return theta1;
}

double logpost_b(vec b, vec h, int T, vec y_T){
  //b = (b0, delta, b2)
  //y = (y0. y1, ..., yT)
  
  double L = 0.0;
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
     
  vec z = y_T.subvec(1, T) - b0 - b1 * y_T.subvec(0, T-1) - b2 * exp(h);
  vec u = exp( -h ) % z;
  
  L -= 0.5 * dot(z, u);
  //# priori b0  
  L += - 0.5 * ((b0 - mu_b0)/s_b0) * (b0 - mu_b0)/s_b0;
  //# priori b1
  L += (a_b1 - 1) * log(1 + b1) + (b_b1 - 1) * log(1 - b1);
  //# priori b2
  L += - 0.5 * ((b2 - mu_b2)/s_b2) * ((b2 - mu_b2)/s_b2);
  //# jacobiano
  L += log(1 - b1 * b1 );
  
  return L;
}

vec glogpost_b(vec b, vec h, mat inv_G, mat dG_b0, mat dG_delta, mat dG_b2, int T, vec y_T){
  //b = (b0, delta, b2)
  // y = (y0. y1, ..., yT)
  
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
  
  vec grad = zeros<vec>(3, 1);
  
  // construindo os vetores u, v e z
  vec z = y_T.subvec(1, T) - b0 - b1 * y_T.subvec(0, T-1) - b2 * exp( h );
  vec v = exp( -h );
  vec u = exp( -h ) % y_T.subvec(0, T-1);
  
  grad[0] += dot(v, z) - (b0 - mu_b0) / (s_b0 * s_b0);
  
  grad[1] += (1 - b1 * b1 ) * ( dot(u, z) ) + (a_b1 - 1) * (1 - b1) - (b_b1 - 1) * (1 + b1); 
  //# jacobiano
  grad[1] -=  2 * b1;
  
  grad[2] += sum( z ) - (b2 - mu_b2) / (s_b2 * s_b2);
  
  // Partial Energy function
  grad[0] += ( -0.5 * trace( inv_G * dG_b0 ) );
  grad[1] += ( -0.5 * trace( inv_G * dG_delta ) );
  grad[2] += ( -0.5 * trace( inv_G * dG_b2 ) );
  
  return - grad;
}


double G_b(vec b, vec h, mat &G, mat &inv_G, mat& dG_b0, mat &dG_delta, mat &dG_b2, int T, vec y_T){
   
  double b1 = tanh( b[1] );
   
  G = zeros<mat>(3, 3);
  inv_G = zeros<mat>(3, 3);
  dG_b0 = zeros<mat>(3, 3);
  dG_delta = zeros<mat>(3, 3);
  dG_b2 = zeros<mat>(3, 3);
	
  G(0, 0) = sum( exp( - h ) ) + 1 / ( s_b0 * s_b0 );
  G(2, 2) = sum( exp( h ) ) + 1 / ( s_b2 * s_b2 );
  G(1, 2) = ( 1 - b1 * b1 ) * sum( y_T.subvec(0, T - 1 ) );
  G(2, 1) = G(1, 2);
    
  for( int j = 0 ; j < T ; j++ ){
      G(1, 1) += ( 1 - b1 * b1 ) * ( 1 - b1 * b1 ) * exp( - h( j ) ) * y_T( j ) * y_T( j );
      G(0, 1) += ( 1 - b1 * b1 ) * exp( - h( j ) ) * y_T( j );
  }

  dG_delta(1, 1) = - 4 * b1 * G(1, 1) - 2 * b1 * ( a_b1 + b_b1 ) * ( 1 - b1 * b1 );
  G(1, 1) += ( a_b1 + b_b1 ) * ( 1 - b1 * b1 );
  G(1, 0) = G(0, 1);
  G(0, 2) = T;
  G(2, 0) = T;

  inv_G = inv_sympd( G );

  dG_delta(0, 1) = - 2 * b1 * G(0, 1);
  dG_delta(1, 0) = dG_delta(0, 1);
  dG_delta(1, 2) = - 2 * b1 * G(1, 2);
  dG_delta(2, 1) = dG_delta(1, 2);

  return 1;
}

vec nuH_b(mat dG_b0, mat dG_delta, mat dG_b2, mat inv_G, vec p){
    vec iGp = inv_G * p;
    vec aux1 = iGp.t() * dG_b0 * iGp;
    vec aux2 = iGp.t() * dG_delta * iGp;
    vec aux3 = iGp.t() * dG_b2 * iGp;
    vec nuH = zeros<vec>(3, 1);
    nuH[0] = aux1[0];
    nuH[1] = aux2[0];
    nuH[2] = aux3[0];
    
    return nuH;
}

vec gradHmom_b(mat inv_G, vec p){
    return inv_G * p;
}

vec rmhmc_b(vec b_cur, vec h, int fixp, int L, double eps, int T, vec y_T , int &acc){
    
    vec theta2p = zeros<vec>(3, 1);
    vec b = zeros<vec>(3, 1);
    int jfix = 1;
    int jgl = 1;
    mat inv_G = zeros<mat>(3, 3);
    mat dG_b0 = zeros<mat>(3, 3);
    mat dG_delta = zeros<mat>(3, 3);
    mat dG_b2 = zeros<mat>(3, 3);
    mat G = zeros<mat>(3, 3);
    
    G_b( b_cur, h, G, inv_G, dG_b0, dG_delta, dG_b2, T, y_T );
    vec pcur = chol( G ).t() * mvrgaussian(3);
   
    vec aux = pcur.t() * inv_G * pcur;
    vec pa = zeros<vec>(3, 1);
    vec pb = zeros<vec>(3, 1);
    vec pc = pcur;
    vec gdpth2(3);
    vec gdmom1(3);
    vec gdmom2(3);
    
    double uMH = 0.0;
    double aMH = 0.0;

    double HMtheta1 = - logpost_b( b_cur, h, T, y_T) + 0.5 * log( det( G )) + 0.5 * aux[0];
    double HMthetaf = 0.0;

    vec gradHMCth2 = zeros<mat>(3, 1);
    vec theta2c = b_cur;
    vec theta2a;
    vec theta2b;
    theta2a = theta2c;
    pa = pc;
    
    for (int jgl = 1 ; jgl < L + 1 ; jgl++ ){

    	gdpth2 = glogpost_b( theta2a, h, inv_G, dG_b0, dG_delta, dG_b2, T, y_T );
    	gradHMCth2 = gdpth2 - 0.5 * nuH_b(dG_b0, dG_delta, dG_b2, inv_G, pa);
    
    	for(int jfix = 1 ; jfix < fixp + 1 ; jfix++ ){
    		pb = pa - 0.5 * eps * gradHMCth2;
     		gradHMCth2 = gdpth2 - 0.5 * nuH_b( dG_b0, dG_delta, dG_b2, inv_G, pb);
    	}
    
    	G_b( theta2a, h, G, inv_G, dG_b0, dG_delta, dG_b2, T, y_T);
    	gdmom1 = gradHmom_b( inv_G, pb );
    	gdmom2 = gradHmom_b( inv_G, pb );
    
    	for(int jfix = 1 ; jfix < fixp + 1 ; jfix++ ){
    		theta2b = theta2a + 0.5 * eps * gdmom1 + 0.5 * eps * gdmom2;
    		G_b( theta2b, h, G, inv_G, dG_b0, dG_delta, dG_b2, T, y_T);
    		gdmom2 = gradHmom_b( inv_G, pb );
    	}
    
    	gdpth2 = glogpost_b( theta2b, h, inv_G, dG_b0, dG_delta, dG_b2, T, y_T );
    	gradHMCth2 = gdpth2 - 0.5 * nuH_b( dG_b0, dG_delta, dG_b2, inv_G, pb);
    	pb = pb - 0.5 * eps * gradHMCth2;
    	pa = pb;
    	theta2a = theta2b;
    }

    aux = pb.t() * inv_G * pb;

    HMthetaf = - logpost_b( theta2b, h, T, y_T ) + 0.5 * log( det( G ) ) + 0.5 * aux[0];
    uMH = randu( distr_param( 0, 1 ) );

    aMH = std::min( 1.0, exp( - HMthetaf + HMtheta1 ) );

    if( uMH < aMH ){
    	b = theta2b;
    	acc++;
    }else{
    	b = b_cur;
    }

    return b;
}


double logpost_h(vec h, vec theta, vec b, int T, vec y_T){
  //h = (h1, ..., hT)
 
  double L = 0.0; 
  double mu = theta[0];
  double phi = tanh( theta[1] );
  double sigma = exp( theta[2] );
  double b0 = b[0];
  double b1 = tanh( b[1] );
  double b2 = b[2];
  
 
  L +=  - 0.5 * sum( h );
 
  for( int j = 1 ; j < T ; j++ ){
    L += - 0.5 / ( sigma * sigma ) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
    L += - 0.5 * exp( -h[j] ) * ( y_T[j+1] - b0 - b1 * y_T[j] - b2 * exp( h[j]) ) * ( y_T[j+1] - b0 - b1 * y_T[j] - b2 * exp(h[j]) );
  }
  
  L += - 0.5 * exp( -h[0] ) * ( y_T[1] - b0 - b1 * y_T[0] - b2 * exp( h[0]) ) * ( y_T[1] - b0 - b1 * y_T[0] - b2 * exp(h[0]) );
  L += - 0.5 * (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma);
  
  return L;
}

vec glogpost_h(vec h, vec theta, vec b, int T, vec y_T){
  //h = (h1, ..., hT)
    
  double mu = theta[0];
  double phi = tanh( theta[1] );
  double sigma = exp(theta[2]);
  double b0 = b[0];
  double b1 = tanh( b[1] );
  double b2 = b[2];
  
  //# construindo o vetor s
  vec r  = zeros<vec>(T, 1);
  vec mu_t = y_T.subvec(1, T) - b0 - b1 * y_T.subvec(0, T - 1) - b2 * exp( h );
  
  vec s = - 0.5 + 0.5 * exp( -h ) % mu_t % mu_t + b2 * mu_t;
  
  //# construindo o vetor r
  r[0] += ( h[0] - phi * h[1] - mu * (1 - phi) ) / (sigma * sigma);
  r[T-1] = 1 / (sigma * sigma) * (h[T-1] - mu - phi * (h[T-2] - mu) );
  
  vec u = h.subvec(1, T - 2);
  vec v = h.subvec(2, T - 1) + h.subvec(0, T - 3);
  
  r.subvec(1, T-2) += (1 + phi * phi ) * u - phi * v - mu * (1 - phi) * (1 - phi);
  r.subvec(1, T-2) /= sigma * sigma ;
  
  return s - r;
}

vec hmc_h(vec h_cur, vec theta, vec b, int L, double eps, int T, vec y_T, int &acc){

    vec vH = zeros<vec>(T , 1);
    vec pcur = mvrgaussian( T );
   
    double Hm1 = 0.0;
    double Hmf = 0.0;
    double uMH = 0.0;
    double aMH = 0.0;
   
    vec auxp = pcur.t() * pcur;
    Hm1 = - logpost_h( h_cur, theta, b, T, y_T) + 0.5 * auxp[0];
    
    vec pa = pcur;
    vec pb = pcur;
    vec vHa = h_cur;
    vec vHb = h_cur;
    
    for (int jlp = 1 ; jlp < L + 1 ; jlp++ ){
    
    	pb = pa + 0.5 * eps * glogpost_h( vHa, theta, b, T, y_T);
    	vHb = vHa + eps * pb;
    	pb = pb + 0.5 * eps * glogpost_h( vHb, theta, b, T, y_T);
    	vHa = vHb;
    	pa = pb;
    	
    }
    
    auxp = pb.t() * pb;
    Hmf = - logpost_h( vHb, theta, b, T, y_T) + 0.5 * auxp[0];

    aMH = std::min( 1.0, exp( - Hmf + Hm1 ) );
    uMH = randu( distr_param( 0, 1 ) );
    
    if( uMH < aMH ){
    	vH = vHb;
    	acc++;
    }else{
        vH = h_cur;
    }
    

    return vH;
}
