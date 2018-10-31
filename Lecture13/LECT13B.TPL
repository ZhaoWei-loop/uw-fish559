DATA_SECTION
 // init in DATA_SECTION means "read it in". These do not carry differentiable information.
 init_int Fyear; // first year
 init_int Lyear; // last year
 init_vector Catch(Fyear,Lyear); // vector indexed from first year to last year
 init_vector CPUE(Fyear,Lyear);
 !! cout << Catch << endl;                                       // Check that the data are read in

PARAMETER_SECTION
// init in the PARAMETER_SECTION means "estimate this parameter"
 init_number logr;
 init_number logq;
 init_number logK;
 number r; // local variables that are available to everywhere, not estimated
 number q;
 number K;
 vector Bio(Fyear,Lyear+1);
 
 objective_function_value objn; // what we're trying to minimize

INITIALIZATION_SECTION 
// one way to do starting values. notice no semicolons in this section.
 logr -0.6
 logq -9
 logK 8.5 

PROCEDURE_SECTION
 int Year;
 dvariable SS;
 
 // Convert from log to normal space
 r = mfexp(logr);
 q = mfexp(logq);
 K = mfexp(logK);

 // Project the model forward
 Bio(Fyear) = K;
 for (Year=Fyear;Year<=Lyear;Year++)
  {
   Bio(Year+1) = Bio(Year) + r*Bio(Year)*(1.0-Bio(Year)/K) - Catch(Year);
   if (Bio(Year+1) < 1) Bio(Year+1) = 1;
  } 
  
 // Compute the sum of squares  
 SS = 0;
 for (Year=Fyear;Year<=Lyear;Year++)
  if (CPUE(Year) > 0)
   SS += square(log(q*Bio(Year)/CPUE(Year)));
 objn = SS;  
 cout << SS << endl;
 

