setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Scripts")
data <- read.table("LectA2.dat", header=TRUE) #data stored in data file
parameters <- list(b0=0, b1=0, logSigma=0) #initialize these params at 0

require(TMB)
compile("LectA2.cpp", flags="-Wno-ignored-attributes") #compile; "who-ignored-attributes" removes warnings so that finding errors is easier
dyn.load(dynlib("LectA2")) #load

################################################################################

model <- MakeADFun(data, parameters, DLL="LectA2",silent=T) #create model object
fit <- nlminb(model$par, model$fn, model$gr) #fit model

best <- model$env$last.par.best
rep <- sdreport(model)

print(best)
print(rep) #maximum gradient vector should be v small (ab partial derivatives being 0); E-05 is pretty good
