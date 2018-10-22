#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(m);
  DATA_IVECTOR(TT);
  DATA_INTEGER(Tmax);
  // matrix where each col is a year, each row is a stock
  DATA_MATRIX(B);
  DATA_MATRIX(R);
  DATA_VECTOR(Phi0);

  // End of data section

  PARAMETER(dummy);
  PARAMETER(mu);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(B0);
  PARAMETER_VECTOR(log_sigR);
  PARAMETER_VECTOR(eps);
  // End of the estimated parameters

  vector<Type> R0(m);
  vector<Type> SigR(m);
  vector<Type> h(m);
  Type tau;
  Type beta;
  Type BH;
  Type obj_fun;
  matrix<Type> log_Rhat(m,Tmax);
  matrix<Type> log_R(m,Tmax);
  //matrix<Type> resids(m,Tmax);
  // End of the temporary variables

  // Transform the parameters
  R0 = B0/Phi0;
  tau = exp(log_tau);
  SigR = exp(log_sigR);
  
  tau = 0;
  beta = 0;
  //obj_fun = 0;
  log_Rhat.setZero();
  log_R.setZero();
  
  for (int k=0;k<m;k++) {
    for (int y=0;y<TT(k);y++) {
      log_R(k,y) = log(R(k,y));
    }
  }
  
  //resids.setZero();
  BH = 0;
  
 for (int k=0;k<m;k++) {
   // Extract beta and define h
   beta=mu+tau*eps(k);
   h(k)=(exp(beta)+Type(0.2))/(Type(1)+exp(beta));
   
   for (int y=0;y<TT(k);y++) {
     
     //Beverton-Holt model
     log_Rhat(k,y)=log(Type(4)*R0(k)*h(k)*B(k,y))-log(Phi0(k)*R0(k)*(Type(1)-h(k)) + 
       (Type(5)*h(k)-Type(1))*B(k,y)); 

     //Likelihood
     BH += pow(log_R(k,y)-log_Rhat(k,y)+(SigR(k)*SigR(k)/Type(2)),Type(2)) /
       (Type(2)*SigR(k)*SigR(k))+log(SigR(k));
     
   }
   //Adjust likelihood for random effect
   BH -= dnorm(eps(k), Type(0), Type(1), true);
   }

  //obj_fun += dummy*dummy;
 
  REPORT(beta);
  REPORT(h);
  REPORT(log_Rhat);
  ADREPORT(h);
  ADREPORT(R0);
  ADREPORT(tau);
  
  //return(obj_fun);
  return(BH);
}
