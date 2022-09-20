#include <TMB.hpp>


template<class Type>
Type posfun(Type x, Type eps, Type &pen) //writing a function; &pen allows that value to be changed outside of fct
{
  if ( x >= eps ){
    return x;
  } else {
    pen += Type(0.01) * pow(x-eps,2);
    return eps/(Type(2.0)-x/eps);
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(C); //C and I both vectors of real numbers
  DATA_VECTOR(I);
  int n = C.size(); //n is length of C

  PARAMETER(logR); //ensuring they cannot be negative
  PARAMETER(logK);
  PARAMETER(logQ);
  PARAMETER(logSigma);
  PARAMETER_VECTOR(FF);
  Type r = exp(logR); //r is real number, assigning exponential of logR
  Type k = exp(logK);
  Type q = exp(logQ);
  Type sigma = exp(logSigma);
  int n1 = 0;
  n1 = n + 1;
  vector<Type> B(n1); //these are all temporary variables; B has length n1
  vector<Type> Ihat(n);
  vector<Type> Chat(n);
  vector<Type> ExpOut(n);
  Type f;
  B(0) = k; //0 is first element of array
  for(int t=0; t<n; t++) //looping over t to n-1
  {
    Type Expl = 1.0/(1.0+exp(-FF(t))); //avoiding "ifs"; expl is logit function of FF, makes it 0 to 1; avoids putting bounds on params; function will never go negative
    B(t+1) = B(t) + r*B(t)*(1-B(t)/k) - Expl*B(t);
    Chat(t) = Expl*B(t);
    ExpOut(t) = Expl;
    Ihat(t) = q*B(t);
  }
  f = -sum(dnorm(log(C), log(Chat), Type(0.05), true));
  f -= sum(dnorm(log(I), log(Ihat), sigma, true)); //f minus equals; equals f minus this thing here
    //take f you had before, substract that quantity
  //f = -sum(dlnorm(C, log(Chat), Type(0.05), true));
  //f -= sum(dlnorm(I, log(Ihat), sigma, true));

  ADREPORT(log(B)); // uncertainty, can print variances of log biomass
  REPORT(f);        
  REPORT(B);        
  REPORT(Chat);        
  REPORT(ExpOut);        
  REPORT(Ihat);     

  return f;
}




// dvariable posfun(const dvariable&x,const double eps,dvariable& pen)
//     {
//       if (x>=eps) {
//         return x;
//       } else {
//         pen+=.01*square(x-eps);
//         return eps/(2-x/eps);
// } }


