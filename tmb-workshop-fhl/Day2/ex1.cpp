#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(model);
  DATA_VECTOR(age);
  DATA_VECTOR(len);
  
  int N = age.size(); // declare and specify an integer n the length of y

  
  PARAMETER(logLinf);
  PARAMETER(logK);
  PARAMETER(logSigma);
  PARAMETER(loga0);
  PARAMETER(loga50);
  PARAMETER(logdelta);
  
  Type Linf = exp(logLinf);
  Type Sigma = exp(logSigma);
  Type K = exp(logK);
  Type a0 = exp(loga0);  
  Type a50 = exp(loga50);
  Type delta = exp(logdelta);
  
  vector<Type> Lpred(N); // create a vector of n real numbers
  vector<Type> resids(N); // create a vector of n real numbers
  vector<Type> ressq(N); // create a vector of n real numbers
  
  if(model == 1) {
  
  Lpred = Linf / (1.0 + exp(-log(19)*(age - a50)/delta));
  resids = len-Lpred;
  ressq = square(resids);
  
  }
  else {
  
  Lpred = Linf * (1.0 - exp(-K * (age - a0)));
  resids = len-Lpred;
  ressq = square(resids);
  
  }
  
  Type neglogL = 0.0; // real number
  
  ADREPORT(Linf);
  ADREPORT(K);
  ADREPORT(Sigma);
  ADREPORT(a0);  
  ADREPORT(a50);  
  ADREPORT(delta);  
  
  //neglogL = (N * logSigma) + sum(square(len-Lpred)) / 2.0*square(Sigma); //0.5 * N * log(2.0 * PI) +
  neglogL = (N * logSigma) + sum(ressq) / (Type(2.0)*pow(Sigma, 2)); //0.5 * N * log(2.0 * PI) +
  //neglogL = -sum(dnorm(len, Lpred, exp(logSigma), true));
  
  return(neglogL);
  //return neglogL;
}
