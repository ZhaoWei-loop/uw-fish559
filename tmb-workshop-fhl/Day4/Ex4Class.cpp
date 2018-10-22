#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nyear);
  DATA_INTEGER(Nclass);
  DATA_VECTOR(Length);
  DATA_VECTOR(Weight);
  DATA_MATRIX(X);
  DATA_VECTOR(S);
  DATA_VECTOR(SurveyS);
  DATA_SCALAR(M);
  DATA_VECTOR(CWObs);
  DATA_MATRIX(CALObs);
  DATA_SCALAR(Neff);
  DATA_VECTOR(BioIndex);
  //DATA_SCALAR(BioSig);
  DATA_INTEGER(Nproj);
  DATA_SCALAR(Fproj);
  DATA_SCALAR(Jac_correction);

  // End of data section

  //PARAMETER(dummy);
  PARAMETER(LogRbar);
  PARAMETER_VECTOR(LogNinit);
  PARAMETER_VECTOR(LogFullF);
  PARAMETER_VECTOR(Eps);

  matrix<Type> N(Nyear+Nproj+1,Nclass);
  matrix<Type> F(Nyear+Nproj,Nclass);
  matrix<Type> Z(Nyear+Nproj,Nclass);
  matrix<Type> CAL(Nyear+Nproj,Nclass);
  vector<Type> CW(Nyear+Nproj);
  vector<Type> BioPred(Nyear+Nproj);

  //Type Penal;
  // Type LikeCatch;
  // Type LikeBio;
  // Type LikeCAL;
  Type obj_fun;
  vector<Type> FullF = exp(LogFullF);
  
  // End of specifications section
  // =============================

  // First set F and Z by size-classs (note that Fproj applies after year Nyear)
  for (int Iyear=0; Iyear<Nyear+Nproj; Iyear++)
  {
   for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
     if (Iyear>=Nyear) 
     {
       F(Iyear,Iclass)=Fproj*S(Iclass);
     } else {
     F(Iyear,Iclass)=FullF(Iyear)*S(Iclass);
     }
     Z(Iyear,Iclass)=M+F(Iyear,Iclass);
    }
  }
  // Set the N matrix
  N.setZero();
  
  // Initialize first year class
  
  for (int Iclass=0;Iclass<Nclass;Iclass++) 
  {   
    //N(Iyear,Iclass)=0;
    N(0,Iclass) = exp(LogNinit(Iclass));
  }
  
  // Fill in the rest of the Catch-at-length and Numbers matricies
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++)
  {
    for (int Iclass=0;Iclass<Nclass;Iclass++)
      {
    // Catch-at-length
    CAL(Iyear,Iclass)=F(Iyear,Iclass)/Z(Iyear,Iclass)*N(Iyear,Iclass)*(Type(1.0)-exp(-Z(Iyear,Iclass)));
      
      for (int IIclass=Iclass;IIclass<Nclass;IIclass++)
      {
        // Numbers-at-length II class is receiving, I is giving
        N(Iyear+1,IIclass)+=X(Iclass,IIclass)*exp(-Z(Iyear,Iclass))*N(Iyear,Iclass);
      }
   }
    // Recruitment (watch fFor the index for Eps - and N)
    N(Iyear+1,0) += exp(LogRbar)*exp(Eps[Iyear]);
  }    

  // Likelihood
  // 0: Catch
  // 1: Biomass
  // 2: CAL (Catch-at-length)
  // 3: Recruitment penalty
  // 4: Jacobian correction for parameters estimated in log space
  
  vector<Type> nll(4);
  nll.setZero();
  
  // Catch Likelihood
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++)
  {
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      CW(Iyear)+=CAL(Iyear,Iclass)*Weight(Iclass);
    }
  }

  nll(0) = -sum(dnorm(log(CWObs), log(CW.head(Nyear)), Type(0.05), true));
  
  // Index Likelihood
  
  // Biomass predictions
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++)
  {
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      BioPred(Iyear)+=SurveyS(Iclass)*N(Iyear,Iclass)*Weight(Iclass);
    }
  }

  // use .sum() because sum() doesn't like vector
  Type logq;
  logq = Type(0);
  logq = (log(BioIndex) - log(BioPred.head(Nyear))).sum() / Nyear;
    
  nll(1) = -sum(dnorm(log(BioIndex), logq + log(BioPred.head(Nyear)), Type(0.2), true));

  // CAL Likelihood
  
  matrix<Type> rho_hat(Nyear,Nclass);
  rho_hat.setZero();
  
  vector<Type> tot_CAL(Nyear);
  tot_CAL.setZero();
  tot_CAL = CAL.topRows(Nyear).rowwise().sum();
  
  for (int Iyear=0;Iyear<Nyear;Iyear++)
  {
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      rho_hat(Iyear,Iclass) = CAL(Iyear,Iclass)/tot_CAL(Iyear);
    }
  }
  
  // Same thing for CALObs
  matrix<Type> rho(Nyear,Nclass);
  rho.setZero();
  
  vector<Type> tot_CALObs(Nyear);
  tot_CALObs.setZero();
  tot_CALObs = CALObs.rowwise().sum();
  
  for (int Iyear=0;Iyear<Nyear;Iyear++)
  {
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      rho(Iyear,Iclass) = CALObs(Iyear,Iclass)/tot_CALObs(Iyear);
    }
  }
 
 for (int Iyear=0;Iyear<Nyear;Iyear++)
 {
   for (int Iclass=0;Iclass<Nclass;Iclass++)
   {
  // nll(2)-=Neff*rho(Iyear,Iclass)*log(rho_hat(Iyear,Iclass))-log(rho(Iyear,Iclass));
  nll(2)-=Neff*CALObs(Iyear,Iclass)*log(rho_hat(Iyear,Iclass))-log(CALObs(Iyear,Iclass));
   }
 }
 
  // Recruitment penalty (include years after Nyear)
  nll(3) = -sum(dnorm(Eps,Type(0),Type(0.6), true));

  obj_fun = nll(0)+nll(1)+nll(2)+nll(3)-Type(3266);//dummy*dummy + 

  // Jacobian correction - check this with STAN example -think it might actually
  // be -par instead of -logpar for the correction
  
  // if (Jac_correction == 1) {
  //  nll(4) -= LogRbar; 
  //   nll(4) -= sum(LogNinit); 
  //   nll(4) -= sum(LogFullF); 
  //   }

  // Stuff to report
  REPORT(N);
  REPORT(F);
  REPORT(BioPred);
  REPORT(obj_fun);
  REPORT(CAL);
  REPORT(CW);
  REPORT(nll);
  REPORT(logq);
  REPORT(rho);
  REPORT(rho_hat);

  return(obj_fun);
}
