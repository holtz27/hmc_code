//###################################################################
//########################## Model ##########################

//function pointer
typedef double ( * num_ptr )(vec, List);
typedef vec    ( * vec_ptr )(vec, List);

// hmc function
mat lf(double eps, int L, vec theta_current, vec p_current, mat inv_M, 
       List param, vec_ptr fun){
  // M = M.i()
  int T = theta_current.n_elem;
  vec theta = theta_current, p = p_current, y = param["y_T"];
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
  double E;
  vec y = param["y_T"];
  
  E = 0.5 * dot(p, inv_M * p);
  E += - fun(theta, param);
  
  return E;
}

List Hmc_one(double eps, int min_L, int max_L, vec theta_init, List param, 
            num_ptr fun1, vec_ptr fun2){
  
  //wall_clock timer;
  //timer.tic();
  
  int acc = 0, dim = theta_init.n_elem; //a = floor(0.1*N)
  double h; 
  
  //inicialindo a cadeia
  
  vec mean(dim, fill::zeros), theta_current = theta_init, p_current;
  mat M(dim, dim, fill::eye), inv_M = M.i(), Proposal(dim, 2);
  
  //(i) gerando p_current e u
  p_current = mvnrnd(mean, M);
    
  //(ii) Leapfrog 
  int L = randi( distr_param( min_L, max_L ) );
  Proposal = lf(eps, L, theta_current, p_current, inv_M, param, fun2);
    
  // (iii) Calculando H's
  h = H(theta_current, p_current, param, inv_M, fun1);
  h -= H(Proposal.col(0), Proposal.col(1), param, inv_M, fun1);
    
  //(iv) Prop aceitação
  if( R::runif(0, 1) < std::min( 1.0, exp(h) ) ){
    theta_current = Proposal.col(0);
    acc += 1;
    }
    
  //Progress
  //if( (i % a) == 0 ) cout << "Progresso em " << ceil(100*i/N)<<" %"<< endl;
  
  
  //double time = timer.toc();
  
  return List::create(Named("theta_current") = theta_current, Named("acc") = acc);
}
// [[Rcpp::export]]
List hmc(double eps, int min_L, int max_L, vec theta_init, List param){
  
  List out = Hmc_one(eps, min_L, max_L, theta_init, param, 
                    &logpost_h, &glogpost_h);
  return out;
  
}
