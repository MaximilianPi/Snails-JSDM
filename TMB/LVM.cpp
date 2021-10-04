// Author: Max
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{ 
  using namespace density;
  DATA_MATRIX(x); // sites*Env
  DATA_MATRIX(y); // sites*species
  int s = y.cols(); // number of species
  int n = x.rows(); // number of sites

  PARAMETER_MATRIX(W); // species:environment response (betas) -> sp*p
  PARAMETER_VECTOR(LF); // species:latent variable response (factor loadings) -> sp*l
  PARAMETER_MATRIX(LV); // latent variables -> n*l
  int l = LV.cols(); // number of latent
  int LV_size = LF.size(); // number of factor loadings
  matrix<Type> fit(n, s); // prediction matrix -> n*sp
  matrix<Type> LF_constrained(l, s); // -> sp*l
  parallel_accumulator<Type> nll(this); // enable parallelization
  Type one=Type(1.0); 
  vector<Type> prob; // predictions on the linear scale
  
  // prior
  for(int j = 0; j < l; j++) {
    nll -= sum( dnorm( vector<Type>(LV.col(j)), Type(0), Type(1), true) );
  }
  
  for(int i = 0; i < LV_size; i++) {
    nll -= ( dnorm( (Type)LF(i), Type(0), Type(1), true) );
  }
  
  for(int i = 0; i < s; i++) {
    nll -= sum( dnorm( vector<Type>(W.col(i)), Type(0), Type(2), true) );
  }
  
  /* construct factor loading matrix with the following constrains:
      - upper triangular = 0 (doesn't work for high l? check Francis Hui's paper again)
      - diagonal in [0, 1]
      - rest in [-1, 1]
  */
  int counter = 0;
  for(int j = 0;j < l; j++) {
    for(int i = 0; i < s; i++) {
      if(j>i) { LF_constrained(i, j) = 0; }
      else if(j == i){ 
        LF_constrained(j,i) = pnorm(LF(counter), Type(0),Type(1)); 
        counter++;
      } else {
        LF_constrained(j,i) = Type(-1.0) + Type(2.0)* pnorm(LF(counter), Type(0),Type(1)); 
        counter++;
      }
    }
  }
  
  // predictions
  fit = x*W + LV*LF_constrained;
  
  // likelihood
  for(int i=0;i<s;i++) {
    prob = vector<Type>( fit.col(i) );
    nll -= sum(  dbinom( vector<Type>(y.col(i)) ,Type(1.0) , vector<Type>( one/(one+exp(-prob)) ) , true)); 
  }
  
  ADREPORT(W);
  
  return nll;
}