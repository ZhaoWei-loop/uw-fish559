#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  int n = y.size(); // declare and specify an integer n the length of y

  PARAMETER(b0);
  PARAMETER(b1);
  PARAMETER(logSigma);
  vector<Type> yfit(n); // create a vector of n real numbers

  REPORT(b0); // report the quantity of b0

  Type neglogL = 0.0; // real number

  yfit = b0 + b1*x;
  neglogL = -sum(dnorm(y, yfit, exp(logSigma), true));

  return neglogL;
}
