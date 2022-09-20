#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x); //x is a vector
  PARAMETER(mu);  //mu is parameter
  PARAMETER(logSigma); //logSigma is parameter
  // this is comment

  Type f;
  f = -sum(dnorm(x,mu,exp(logSigma), true)); //f is negative summation of normal distribution w data x, mean mu,
                                             //and variance logSigma
                                             //if log is true, probabilities p are given as log(p)
  REPORT(mu);
  return f; //returns f but also reports mu
}

