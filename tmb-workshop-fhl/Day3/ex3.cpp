#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA SECTION
  DATA_INTEGER(nyears);
  DATA_INTEGER(njourn);
  DATA_INTEGER(model);
  DATA_MATRIX(data);

  // PARAMETER SECTION
  PARAMETER_VECTOR(b0);
  PARAMETER_VECTOR(log_b1);
  PARAMETER_VECTOR(amp);
  PARAMETER_VECTOR(phase);
  PARAMETER_VECTOR(log_period);
  
  vector<Type> period = exp(log_period);
  vector<Type> b1 = exp(log_b1);

  // MODEL OBJECTS
  Type obj_fun = 0;
  matrix<Type> cite_hat(nyears, njourn);
  
  //if( prey.cols() != consum.cols()){
  //  error("Column size of prey and consumption do not match.");
  //}
  
  // FIT MODELS
  for(int i = 0; i < nyears; i++){
    for(int j = 0; j < njourn; j++){
      
      // Linear model
      switch( model ){
      case 1 : // Holling type I
        cite_hat(i,j) = b0(j)+b1(j)*data(i,0)+amp(j)*sin(Type(2)*PI*(data(i,0)-phase(j))/period(j));
        break;
        
      case 2 : //Exponential model
        cite_hat(i,j) = b0(j)*(1-b1(j)*exp(-data(i,0)))+amp(j)*sin(Type(2)*PI*(data(i,0)-phase(j))/period(j)); 
      //  break;
        
      }
      
      // Fit objective function
      //obj_fun += -log(dpois(data(i,j+1), cite_hat(i,j)));
      obj_fun -= dpois(data(i,j+1), cite_hat(i,j), true);
    }
  }
  
  REPORT(cite_hat);
  REPORT(period);
  REPORT(b1);
  
  return obj_fun;
}
