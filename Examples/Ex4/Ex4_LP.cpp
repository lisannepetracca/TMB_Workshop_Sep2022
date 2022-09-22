#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nyear); //there are 25 years
  DATA_INTEGER(Nclass); //there are 5 age classes
  DATA_VECTOR(Length); //there are 5 lengths
  DATA_VECTOR(Weight); //there are 5 weights
  DATA_MATRIX(X);//X is 5x5 growth transition matrix
  DATA_VECTOR(S);//S is selectivity by 5 age classes
  DATA_VECTOR(SurveyS); //survey selectivity for animals in each size class
  DATA_SCALAR(M); //instantaneous rate of natural mortality
  DATA_VECTOR(CWObs); //observed catch by weight
  DATA_MATRIX(CALObs); //observed proportion of catch by age class
  DATA_SCALAR(Neff); //effective sample size (set to 100)
  DATA_VECTOR(BioIndex); //some form of biomass in vector form
  DATA_SCALAR(BioSig); //another form of biomass
  DATA_INTEGER(Nproj); //projected N
  DATA_SCALAR(Fproj); //projected fishing mortality rate

  // End of data section

  PARAMETER(dummy);
  PARAMETER(LogRbar);
  PARAMETER_VECTOR(LogNinit);
  PARAMETER_VECTOR(LogFullF);
  PARAMETER_VECTOR(Eps);

  matrix<Type> N(Nyear+Nproj+1,Nclass);
  matrix<Type> F(Nyear+Nproj,Nclass);
  matrix<Type> Z(Nyear+Nproj,Nclass);
  matrix<Type> CAL(Nyear+Nproj,Nclass);
  vector<Type> CW(Nyear+Nproj);
  vector<Type> BioPred(Nyear+Nproj);

  Type CALtot;

  Type Penal;
  Type LikeCatch;
  Type LikeBio;
  Type LikeCAL;
  Type obj_fun;

  // End of specifications section
  // =============================

  // First set F and Z by size-classs (note that Fproj applies after year Nyear)
  for (int Iyear=0; Iyear<Nyear+Nproj; Iyear++) //looping over years
    for (int Iclass=0;Iclass<Nclass;Iclass++) //looping over age classes
    {
      if (Iyear < Nyear) //for years in data period
        F(Iyear,Iclass) = exp(LogFullF(Iyear))*S(Iclass); //here is the mortality from fishing
      else //for years in the projection period
        F(Iyear,Iclass) = Fproj*S(Iclass); //here is the mortality from fishing
      Z(Iyear,Iclass) = M + F(Iyear,Iclass); //here is natural and fishing mortality by class and year
    }
    
    // Now set the N matrix
    for (int Iclass=0;Iclass<Nclass;Iclass++) N(0,Iclass) = exp(LogNinit(Iclass)); //setting N to initial
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++) //running thru all years
  {
    // Catch-at-length
    CALtot = 0; CW(Iyear) = 0;
    for (int Iclass=0;Iclass<Nclass;Iclass++) //running over age classes
    {
      CAL(Iyear,Iclass) = F(Iyear,Iclass)/Z(Iyear,Iclass)*N(Iyear,Iclass)*(1.0-exp(-Z(Iyear,Iclass))); //catch at length
      CALtot += CAL(Iyear,Iclass); //adding these sequentially
      CW(Iyear) += Weight(Iclass)*CAL(Iyear,Iclass); //getting catch by weight iteratively for a given year
    }
    CALtot = (CAL.row(Iyear)).sum(); //summing across rows to get total catch at length
    for (int Iclass=0;Iclass<Nclass;Iclass++) CAL(Iyear,Iclass) /= CALtot; //dividing all these by CALtot to get proportions
    
    // Numbers-at-age
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      N(Iyear+1,Iclass) = 0;
      for (int Jclass=0;Jclass<Nclass;Jclass++)
        N(Iyear+1,Iclass) += N(Iyear,Jclass)*exp(-Z(Iyear,Jclass))*X(Jclass,Iclass); //this is main N equation
    }
    
    // Recruitment (watch for the index for Eps - and N)
    N(Iyear+1,0) += exp(LogRbar)*exp(Eps[Iyear]);
  }
  
  
  // Catch Likelihood
  Type SS = 0; //stands for sum of squares
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++) //running thru all years
    SS += square(log(CWObs(Iyear)) - log(CW(Iyear)));
    LikeCatch = SS / (2 * 0.05 * 0.05);
    
    // Biomass predictions
    for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++){ //running thru all years
      BioPred(Iyear) = 0;
      for (int Iclass=0;Iclass<Nclass;Iclass++) {
        BioPred(Iyear) += SurveyS(IClass) * Weight(Iclass) * N(Iyear,IClass); //something ab year starting at +1
    }}
  
  // Index Likelihood
  
  Type Top = 0; Type Bot = 0; Type q; //this is some weird q stuff; bottom is adding 1 for each year; q equation at bottom
  for (int Iyear=0; Iyear<Nyear; Iyear++)
  { Top += log(BioIndex(Iyear)/BioPred(Iyear)); Bot += 1.0; }
  q = exp(Top/Bot);
  
  SS = 0;
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++) //running thru all years
    SS += square(log(BioIndex(Iyear)) - log(BioPred(Iyear)*q)); 
    LikeBio = SS / (2 * BioSig * BioSig); //BioSig represents the SE on log of index of abundance 
                                         //why not set at 0.2?
  //THIS IS WHERE I STOPPED, AND HONESTLY I HAD DONE PRETTY WELL UP UNTIL THIS POINT!!!
  
  // CAL Likelihood
  LikeCAL = 0;
  for (int Iyear=0; Iyear<Nyear; Iyear++)
    for (int Iclass=0;Iclass<Nclass;Iclass++)
      if (CALObs(Iyear,Iclass) > 0) //if observed CAL greater than 0, calculate the likelihood
        LikeCAL -= Neff*CALObs(Iyear,Iclass)*log(CAL(Iyear,Iclass)/CALObs(Iyear,Iclass)); //equation 7 in assignment

  // Recruitment penalty (include years after Nyear)
  Penal = 0;

  obj_fun = dummy*dummy + LikeCatch+LikeBio+LikeCAL+Penal;

  // Stuff to report
  REPORT(N);
  REPORT(LikeCatch);
  REPORT(LikeBio);
  REPORT(LikeCAL);
  REPORT(Penal);
  REPORT(BioPred);
  REPORT(obj_fun);

  return(obj_fun);
}
