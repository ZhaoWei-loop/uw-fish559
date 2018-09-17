#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nyear)
  DATA_INTEGER(Nage)
  int Nlast = Nage-1;
  DATA_SCALAR(M)
  DATA_VECTOR(Wght)
  DATA_VECTOR(CatchW)
  DATA_VECTOR(CatchN)
  DATA_VECTOR(RecIndex)
  DATA_VECTOR(Effort)
  DATA_INTEGER(RandomR)
  DATA_INTEGER(RandomQ)
  DATA_INTEGER(UseRecIndex)
  // End of data section

  PARAMETER(dummy);
  PARAMETER(F0);
  PARAMETER(LogRec1983);
  Type Rec1983 = exp(LogRec1983);
  PARAMETER_VECTOR(LogRec);
  PARAMETER(logq);
  Type q = exp(logq);
  PARAMETER(LogSigCatch);
  Type SigCatch = exp(LogSigCatch);
  PARAMETER(LogSigmaR);
  Type SigmaR = exp(LogSigmaR);
  PARAMETER(LogSigRecIndex);
  Type SigRecIndex = exp(LogSigRecIndex);

  PARAMETER_VECTOR(logqdev);
  PARAMETER(LogSigmaQ);
  Type SigmaQ = exp(LogSigmaQ);
  PARAMETER(LogitRho);
  Type Rho = -0.999 + Type(2)*0.999/(1+exp(LogitRho));

  // End of the parameter section

  matrix<Type> N(Nyear+1,Nage);
  vector<Type> F(Nyear);
  vector<Type> Z(Nyear);
  vector<Type> Exploit(Nyear);

  vector<Type> Rec(Nyear);
  vector<Type> Bio(Nyear);
  vector<Type> Ntot(Nyear);
  vector<Type> CpredN(Nyear);
  vector<Type> CpredW(Nyear);

  vector<Type> qdev(Nyear);
  vector<Type> qv(Nyear);
  qdev(0) = 0;
  for (int Year=1;Year<Nyear;Year++)
   qdev(Year) = Rho*qdev(Year-1)+sqrt(1.0-Rho*Rho)*logqdev(Year-1);

  // Compute F, Z and exploitation rate
  for (int Year=0;Year<Nyear;Year++)
   {
	qv(Year) = q*exp(qdev(Year));
	F(Year) = qv(Year)*Effort(Year);
	Z(Year) = M + F(Year);
	Exploit(Year) = F(Year)/Z(Year)*(1.0 - exp(-Z(Year)));
   }

  Type Like1;
  Type Like2;
  Type Like3;
  Type obj_fun;

  // End of specifications section
  // =============================

  // Set up in the initial N vector
  N(0,0) = Rec1983;
  for (int Age=1;Age<=Nlast;Age++) N(0,Age) = N(0,Age-1)*exp(-F0-M);
  N(0,Nlast) = N(0,Nlast)/(1.0-exp(-F0-M));

  // Project forward
  for (int Year=0;Year<Nyear;Year++)
   {
    // Basic dynamics
    for (int Age=1;Age<=Nlast;Age++) N(Year+1,Age) = N(Year,Age-1)*exp(-Z(Year));
    N(Year+1,Nlast) = N(Year+1,Nlast) + N(Year,Nlast)*exp(-Z(Year));

    // Totals
    Bio(Year) = 0; Ntot(Year) = 0;
    for (int Age=0;Age<=Nlast;Age++)
     {
	  Bio(Year) +=  Wght(Age)*N(Year,Age);
	  Ntot(Year) += N(Year,Age);
	  Rec(Year) = N(Year,0);
	 }
	CpredN(Year) = Exploit(Year)*Ntot(Year);
	CpredW(Year) = Exploit(Year)*Bio(Year);

    // Recruitment
    if (Year!=Nyear-1)
     N(Year+1,0) = exp(LogRec(Year));
    else
     N(Year+1,0) = Type(1000);
   }
  //std::cout << N << "\n";
  //std::cout << Ntot << "\n";
  //std::cout << Bio << "\n";
  //std::cout << Wght << "\n";

  // Find the ML estimate of q
  Type MeanRecQ = 0;
  for (int Year=0;Year<Nyear;Year++) MeanRecQ += log(Rec(Year)/RecIndex(Year));
  MeanRecQ /= float(Nyear);

  Like1 = 0; Like2 = 0; Like3 = 0;
  Type Residual;
  for (int Year=0;Year<Nyear;Year++)
   {
    Residual = CatchN(Year) - CpredN(Year);
    Like1 += log(SigCatch) + 0.5*square(Residual)/square(SigCatch);
    Residual = CatchW(Year) - CpredW(Year);
    Like2 += log(SigCatch) + 0.5*square(Residual)/square(SigCatch);
    Residual = log(Rec(Year)/RecIndex(Year)) - MeanRecQ;
    Like3 += log(SigRecIndex) + 0.5*square(Residual)/square(SigRecIndex);
   }

  // Find the mean of the rec_devs and the hence the probability of the rec-devs
  Type MeanLogRec = 0;
  for (int Year = 0; Year<Nyear-1;Year++) MeanLogRec += LogRec(Year);
  MeanLogRec /= float(Nyear-1);
  Type ProbRecDev = 0;
  for (int Year = 0; Year<(Nyear-1);Year++)
   ProbRecDev += log(SigmaR) + 0.5*square(LogRec(Year)-MeanLogRec)/square(SigmaR);

  // Only used if needed
  if (RandomR==0) ProbRecDev = 0;

  Type ProbQDev = 0;
  for (int Year = 0; Year<(Nyear-1);Year++)
   ProbQDev += log(SigmaQ) + 0.5*square(logqdev(Year))/square(SigmaQ);

  // Only used if needed
  if (RandomQ==0) ProbQDev = 0;

  // Should the recruitment index be used
  if (UseRecIndex==0) Like3 = 0;

  // Objective function
  obj_fun = dummy*dummy+Like1+Like2+Like3+ProbRecDev+ProbQDev;

  REPORT(N);
  REPORT(Ntot);
  REPORT(CpredN);
  REPORT(CpredW);
  REPORT(Like1);
  REPORT(Like2);
  ADREPORT(qv);
  ADREPORT(Bio);
  ADREPORT(Rec);
  ADREPORT(SigCatch);
  ADREPORT(SigmaR);
  ADREPORT(SigmaQ);
  ADREPORT(SigRecIndex);
  ADREPORT(q);

  return(obj_fun);
}
