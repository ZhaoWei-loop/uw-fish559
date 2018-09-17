#include <TMB.hpp>
 
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  PARAMETER(mu);
  PARAMETER(logSigma);
  
  Type f; // template variable, identifies f as a real number (like 'float' in C++)
  f = -sum(dnorm(x,mu,exp(logSigma), true)); // true moves this to log space
  return f;
}
