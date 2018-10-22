#include <TMB.hpp>


template <class Type> Type square(Type x){return x*x;}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Predator);
  DATA_MATRIX(Prey);
  DATA_MATRIX(Consump);
  DATA_INTEGER(Ndata);
  DATA_INTEGER(Model_type);

  PARAMETER_VECTOR(logalpha);
  PARAMETER_VECTOR(logbeta);
  PARAMETER_VECTOR(gamma);
  PARAMETER_VECTOR(loglambda);

  vector<Type> alpha = exp(logalpha);
  vector<Type> beta = exp(logbeta);
  vector<Type> lambda = exp(loglambda);

  matrix<Type> Pred(Ndata, 3);
  //vector<Type> PredY(3);
  
  //Type dummy = 2;

  Type neglogL = 0.0;

  if (Model_type == 1)
   {
    
    for(int obs = 0; obs < Ndata; obs++) 
    {
      for(int spp = 0; spp < 3; spp++)
      {
      Pred(obs, spp) = alpha(spp) * Predator(obs);
      }}
   }
	//for (II=0;II<=3;II++)
	// PredY(II-1) = Linf/(1+exp(-log(19)*(float(II)-a50)/Delta));
  // }
   
  if (Model_type == 2)
   {
    for(int obs = 0; obs < Ndata; obs++)
    {
      for(int spp = 0; spp < 3; spp++)
      {
    Pred(obs, spp) = (alpha(spp) * Predator(obs))/(1 + beta(spp) * Prey(obs, spp));
      }}
	//for (II=1;II<=20;II++)
	// PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0)));
   }
  
  if (Model_type == 3)
  {
    for(int obs = 0; obs < Ndata; obs++)
    {
      for(int spp = 0; spp < 3; spp++)
      {
    Pred(obs, spp) = (alpha(spp) * Predator(obs) * pow(Prey(obs, spp), gamma(spp) - 1)) / (1 + pow(beta(spp) * Prey(obs, spp), gamma(spp)));
      }}
  //for (II=1;II<=20;II++)
    // PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0)));
  }
  
  if (Model_type == 4)
  {
    for(int obs = 0; obs < Ndata; obs++)
    {
      for(int spp = 0; spp < 3; spp++)
      {
    Pred(obs, spp) = (alpha(spp) * Predator(obs)) / (1 + beta(spp) * Prey(obs, spp) + lambda(spp) * Predator(obs));
      }}
    //for (II=1;II<=20;II++)
    // PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0)));
  }
  
 for(int obs = 0; obs < Ndata; obs++) 
 {
  for(int spp = 0; spp < 3; spp++)
  {
    //neglogL += -sum(dnorm(log(Consump(obs,spp)), log(Pred(obs,spp)), 1, true));
    neglogL += pow(log(Consump(obs,spp)) - log(Pred(obs,spp)),2);
  //std::cout << alpha << " " << neglogL  << "\n";
  }}
 
  REPORT(Pred);
  ADREPORT(alpha);
  //REPORT(dummy);
  return neglogL;
}
