#include <TMB.hpp>


//ok so this example was fitting two types of logistic growth curves to data

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Age); //age and length are vector data
  DATA_VECTOR(Length);
  DATA_INTEGER(Ndata) //Ndata is integer for number of data rows
  DATA_INTEGER(Model_type) //there are two model types

  PARAMETER(LogLinf); //the logs mean they are constrained to be positive
  PARAMETER(Loga50);
  PARAMETER(LogDelta);
  PARAMETER(LogKappa);
  PARAMETER(t0); //this has something to do with time 0 lol
  PARAMETER(LogSigma);

  Type Linf = exp(LogLinf); //bringing back to non-log using Type
  Type a50 = exp(Loga50);
  Type Delta = exp(LogDelta);
  Type Kappa = exp(LogKappa);
  Type Sigma = exp(LogSigma);


  vector<Type> Pred(Ndata); //ok so pred is length of Ndata
  vector<Type> PredY(20);   //PredY only has length of 20
  int II; //II is an integer

  Type neglogL = 0.0; //neglogL starts at 0

  if (Model_type == 1)
   {
	Pred = Linf/(1+exp(-log(19)*(Age-a50)/Delta)); //predictions for which model type is 1
	for (II=1;II<=20;II++) //loops from 1 to 20; note that indices start at 0 so have to subtract 1 below
	 PredY(II-1) = Linf/(1+exp(-log(19)*(float(II)-a50)/Delta)); //here age is from 1:20 it seems
   }
  if (Model_type == 2)
   {
    Pred = Linf*(1-exp(-Kappa*(Age-t0)));
	for (II=1;II<=20;II++)
	 PredY(II-1) = Linf*(1-exp(-Kappa*(float(II)-t0))); //a note that this is predicting for ages 1:20
	                                                    //hence II goes from 1 to 20
	                                                    //the regular preds are not from evenly ordered ages like that
   }
  neglogL = -sum(dnorm(Length, Pred, exp(LogSigma), true)); //negloglik is -1 * sum of normal dist w data length,
                                                            //mean pred, and sd sigma
  //std::cout << Linf << " " << a50 << " " << Delta << " " << exp(LogSigma) << " " << neglogL  << "\n";

  REPORT(PredY);
  REPORT(Pred); //LP added this
  ADREPORT(Linf); //saves standard errors for these things
  ADREPORT(a50);
  ADREPORT(Delta);
  ADREPORT(Kappa);
  ADREPORT(Sigma);


  return neglogL;
}
