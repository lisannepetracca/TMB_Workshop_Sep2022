#include <TMB.hpp>

//this was me doing the first part of the exercise, linear model only
//predicting consumption based on predator densities

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Consump); //consumption is 100 x 3 matrix [1 col per prey]
  DATA_VECTOR(Predator); //predator is vector
  DATA_INTEGER(Ndata);  //Ndata is 100
  DATA_INTEGER(Nprey);  //Nprey is 3
    
  PARAMETER_VECTOR(alpha); //we are estimating slope alpha, and there are three values
  
  matrix<Type> Pred(Ndata, Nprey); //here are a bunch of matrices; first index is row, then col
  matrix<Type> Pred_log(Ndata, Nprey);
  matrix<Type> Data_log(Ndata, Nprey);
  matrix<Type> Diff(Ndata, Nprey);

  int i; //have two loops here
  int j;
  
  Type sumsq = 0.0; //start sumsq at 0
  
    for (i=0;i<Nprey;i++) {
      Pred.col(i) = alpha(i) * Predator; //this is saying the ith column of pred
          for(j=0;j<Ndata;j++){
        Pred_log(i,j) = log(Pred(i,j));
        Data_log(i,j) = log(Consump(i,j));
        Diff(i,j) = Data_log(i,j) - Pred_log(i,j);
        sumsq += pow(Diff(i,j),2);
      }
    }
  
  return sumsq;
}

