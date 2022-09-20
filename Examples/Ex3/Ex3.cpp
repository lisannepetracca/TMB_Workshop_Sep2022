#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Ndata); //ok so this is our number of data rows
  DATA_VECTOR(Google); //this is our data with google citations
  vector<Type> Pred(Ndata); //our predictions have length Ndata
  int i; //i is an integer
  PARAMETER(logLambda); //forcing lambda to be positive
  Type lambda = exp(logLambda);
  
  Type neglogL = 0.0; //starting the negative log lik at 0
  
  neglogL = -sum(dpois(Google, lambda, true)); //equal to negative 1 times the sum
                                               //here we are adding up values of poisson distribution w data
                                               //Google and estimated lambda
  
  REPORT(lambda); //we want lambda reported

  return neglogL; //and neglogL returned
}
