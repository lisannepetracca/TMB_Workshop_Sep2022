#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(m); //number of stocks
  DATA_IVECTOR(TT); //number of years
  DATA_INTEGER(Tmax); //max number of years
  DATA_MATRIX(B); //biomass matrix
  DATA_MATRIX(R); //recruitment matrix
  DATA_VECTOR(Phi0); //spawner biomass per recruit

  // End of data section

  PARAMETER(dummy); 
  PARAMETER(mu); //mean of random effect
  PARAMETER(log_tau); //log of sd of random effect
  PARAMETER_VECTOR(B0); //initial steepness; bounded between 0.1 and 1000
  PARAMETER_VECTOR(log_sigR); //i think this is log of sd of F
  PARAMETER_VECTOR(eta); //no bounds, part of steepness equation 
  // End of the estimated parameters

  vector<Type> R0(m); //stock-specific unfished recruitment; is B0/Phi0
  vector<Type> SigR(m); //this is sd on recruitment
  vector<Type> h(m); //steepness
  vector<Type> beta; //this is random effect of steepness
  matrix<Type> BH(m,Tmax); //this represents the first part of R equation; aka F*Bit 
  matrix<Type> eps(m,Tmax); //i created this to represent the epsilon in the main Rit equation
  
  Type tau; //sd of random effect
  Type obj_fun;
  Type NegLogLike = 0;
  
  // End of the temporary variables

  // Transform the parameters
  R0 = B0/Phi0;
  tau = exp(log_tau);
  SigR = exp(log_sigR);
  obj_fun = 0;
  
  // Extract beta and define h
  for (int i=0;i<m;i++) {
    beta[i] = log((h[i]-0.2)/(1-h[i]));
   for (int j=0;j<TT;j++) {
    BH[i,j] = (4 * R0[i] * h[i] * B[i,j]) / (Phi0[i] * R0[i] * (1-h[i]) + (((5 * h[i]) - 1) * B[i,j]));
    R[i,j] = BH[i,j] * exp(eps[i,j] - (SigR[i]*SigR[i]/2));
    eps[i,j] ~ dnorm(0,SigR);
    B[i,j] ~ dnorm(mu,tau);
   }}

  NegLogLike += (log(R) - log(BH) + ((SigR * SigR)/2)) * (log(R) - log(BH) + ((SigR * SigR)/2)) / ((2 * SigR * SigR) + log_sigR)
  
  obj_fun += dummy*dummy;
  ADREPORT(h);
  ADREPORT(R0);
  ADREPORT(tau);

  return(obj_fun);
}
