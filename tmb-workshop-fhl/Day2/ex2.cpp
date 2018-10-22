#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Predator);
  DATA_VECTOR(Prey);
  DATA_VECTOR(Consump);
  DATA_INTEGER(Ndata);
  DATA_INTEGER(Model_type);

  PARAMETER(logalpha);
  PARAMETER(logbeta);
  PARAMETER(gamma);
  PARAMETER(loglambda);

  Type alpha = exp(logalpha);
  Type beta = exp(logbeta);
  Type lambda = exp(loglambda);

  vector<Type> Pred(Ndata);
  vector<Type> PredY(3);
  
  //Type dummy = 2;

  Type neglogL = 0.0;

  if (Model_type == 1)
   {
	Pred = alpha * Predator;
    
	//for (II=0;II<=3;II++)
	// PredY(II-1) = Linf/(1+exp(-log(19)*(float(II)-a50)/Delta));
   }
   
  if (Model_type == 2)
   {
    Pred = (alpha * Predator)/(1 + beta * Prey);
	//for (II=1;II<=20;II++)
	// PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0)));
   }
  
  if (Model_type == 3)
  {
    Pred = (alpha * Predator * pow(Prey, gamma - 1)) / (1 + pow(beta * Prey, gamma));
    //for (II=1;II<=20;II++)
    // PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0)));
  }
  
  if (Model_type == 4)
  {
    Pred = (alpha * Predator) / (1 + beta * Prey + lambda * Predator);
    //for (II=1;II<=20;II++)
    // PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0)));
  }
  
  neglogL = -sum(dnorm(log(Consump), log(Pred), 1, true));
  std::cout << alpha << " " << neglogL  << "\n";

  REPORT(Pred);
  //REPORT(dummy);
  return neglogL;
}
