#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Length);
  DATA_VECTOR(Age);
  int n = Length.size();
  PARAMETER(logLinf);
  //PARAMETER(logKappa);
  //PARAMETER(a0);
  PARAMETER(a50);
  PARAMETER(logSigma);
  //Type Kappa = exp(logKappa);
  Type Sigma = exp(logSigma);
  Type Linf = exp(logLinf);
  Type neglogL2 = 0.0;
  vector<Type> Pred2(n);
  //Pred = Linf*(1.0-exp(-Kappa*(Age-a0)));
  //neglogL = -sum(dnorm(Length, Pred, Sigma, true));
  //return neglogL;
  
  Pred2 = Linf*(1/(1+exp(-log(19)*((Age-a50)/((0.95*Linf)-a50)))));
  neglogL2 = -sum(dnorm(Length, Pred2, Sigma, true));
  return neglogL2;
  
} 
