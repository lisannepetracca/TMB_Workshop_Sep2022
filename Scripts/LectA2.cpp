#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  int n = y.size(); //n is integer with length y

  PARAMETER(b0);
  PARAMETER(b1);
  PARAMETER(logSigma);
  vector<Type> yfit(n); //yfit is a vector of length n

  REPORT(b0); //we are reporting what i believe is biomass

  Type neglogL = 0.0; //starting neglogL at 0

  yfit = b0 + b1*x; //yfit as linear model
  neglogL = -sum(dnorm(y, yfit, exp(logSigma), true)); //neglogL is -1* sum of all those data points
                                                       //fit to a normal distribution w mean yfit and 
                                                       //sd of logSigma
                                                       //something about being on log scale constraining values to be positive

  return neglogL;
  
}


  