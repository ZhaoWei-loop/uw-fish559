#include <TMB.hpp>

//template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(MM);
  DATA_VECTOR(FF);
  DATA_INTEGER(model);
  int Nyear = MM.size();

  PARAMETER(p);
  PARAMETER(dummy);
  PARAMETER(tau);

  Type logit_p = log(p/(Type(1.0)-p));

  
  vector<Type> pfit(Nyear);
  vector<Type> sigma_y(Nyear);
  //vector<Type> log_sigma_y(Nyear);
  //sigma_y.setZero();
  Type neglogL = 0;
  
  
  for (int i=0;i<Nyear;i++) {
    
    // prediction
    pfit(i) = MM(i)/(MM(i)+FF(i));    
    sigma_y(i) = sqrt(pow(tau,Type(2.0)) + pfit(i)*(Type(1.0)-pfit(i))/(MM(i)+FF(i)));

    switch(model) {
    case 1 : // Model 1 (binomial)
      neglogL -= dbinom_robust(MM(i), MM(i)+FF(i), logit_p, true);
      break;
      
    case 2 : // Model 2 (binomial with overdispersion)
      neglogL += log(sigma_y(i)) + Type(1.0)/(Type(2.0)*pow(sigma_y(i),Type(2.0)))*(pow(pfit(i) - p, Type(2.0)));
      break;
      
    //case 3 : // Model 3 (negative binomial)  
    }
  }
  
  //REPORT(log_sigma_y);
  REPORT(sigma_y);
  REPORT(pfit);
  return neglogL;
}
