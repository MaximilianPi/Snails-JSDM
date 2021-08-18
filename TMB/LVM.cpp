// Author: Max
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{ 
  using namespace density;
  DATA_MATRIX(x);
  DATA_MATRIX(y);
  int s = y.cols();
  int n = x.rows();

  PARAMETER_MATRIX(W);
  PARAMETER_VECTOR(LF);
  PARAMETER_MATRIX(LV);
  int l = LV.cols();
  int LV_size = LF.size();
  matrix<Type> fit(n, s);
  matrix<Type> LF_constrained(l, s);
  parallel_accumulator<Type> nll(this);
  Type one=Type(1.0);
  vector<Type> prob;
  
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
  
  REPORT(W);
  
  return nll;
}