#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x);
  DATA_MATRIX(y);
  int s = y.cols();
  int n = x.rows();
  
  PARAMETER_MATRIX(W);
  PARAMETER_MATRIX(LF);
  PARAMETER_MATRIX(LV);
  matrix<Type> fit(n, s);
  parallel_accumulator<Type> nll(this);
  Type one=Type(1.0);
  vector<Type> prob;
  int l = LV.cols();
  
  // prior
  for(int j = 0; j < l; j++) {
    nll -= sum( dnorm( vector<Type>(LV.col(j)), Type(0), Type(1), true) );
  }
  
  fit = x*W + LV*LF;
  
  for(int i=0;i<s;i++) {
    prob = vector<Type>( fit.col(i) );
    nll -= sum(  dbinom( vector<Type>(y.col(i)) ,Type(1.0) , vector<Type>( one/(one+exp(-prob)) ) , true)); 
  }
  
  REPORT(W);
  
  return nll;
}