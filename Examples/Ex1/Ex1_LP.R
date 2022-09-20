setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Examples/Ex1")
getwd()
data <- read.table("ex1.dat", skip=2, header=T) #data stored in data file
head(data)

require(TMB)
compile("Ex1.cpp", flags="-Wno-ignored-attributes") #compile; "who-ignored-attributes" removes warnings so that finding errors is easier
dyn.load(dynlib("Ex1")) #load

################################################################################
##VERSION 2 OF MODEL

parameters <- list(logLinf=log(20), logKappa=log(4), a50=10, logSigma=log(2)) #a0=0
model <- MakeADFun(data, parameters, DLL="Ex1",silent=T) #create model object
fit <- nlminb(model$par, model$fn, model$gr) #fit model

best <- model$env$last.par.best
rep <- sdreport(model)

# the negative log-likelihood 
minusLL <- function(logLinf,logK,a0,logSigma,Ages,Lengths)
{
  # Extract the parameters
  Linf <- exp(logLinf); Kappa <- exp(logK); Sigma <- exp(logSigma)
  
  # make the model predictions
  Pred <- Linf*(1.0-exp(-Kappa*(Ages-a0)))
  
  # Compute the negative log-likelihood
  NegLogL <- -1*sum(dnorm(Lengths,Pred,Sigma,TRUE))
  
  return(NegLogL)
  
} 