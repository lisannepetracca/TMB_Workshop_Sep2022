setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Scripts")




#==========================================================================================

require(TMB)
compile("LectC3.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("LectC3"))

################################################################################

Nyear <- scan("LectC3.dat",skip=1,n=1,quiet=T)
Nage <- scan("LectC3.dat",skip=3,n=1,quiet=T)+1
M <- scan("LectC3.dat",skip=5,n=1,quiet=T)
Wght <- scan("LectC3.dat",skip=7,n=Nage,quiet=T)
SigCatch <- scan("LectC3.dat",skip=9,n=1,quiet=T)
SigCPUE <-scan("LectC3.dat",skip=11,n=1,quiet=T)  
Omega <-scan("LectC3.dat",skip=13,n=1,quiet=T)
Catch <- matrix(scan("LectC3.dat",skip=15,n=Nyear*3,quiet=T),ncol=3,byrow=T)[,2]
CPUE <- matrix(scan("LectC3.dat",skip=15,n=Nyear*3,quiet=T),ncol=3,byrow=T)[,3]
Propn <- matrix(scan("LectC3.dat",skip=37,n=Nyear*(Nage+1),quiet=T),ncol=Nage+1,byrow=T)[,-1]

################################################################################

LogN <- scan("LectC3.pin",skip=2,n=Nyear+Nage-1,quiet=T)
Sel50 <- scan("LectC3.pin",skip=4,n=1,quiet=T)
Sel95 <- scan("LectC3.pin",skip=6,n=1,quiet=T)
LogFish <- scan("LectC3.pin",skip=8,n=Nyear,quiet=T)
logq <- scan("LectC3.pin",skip=10,n=1,quiet=T)

################################################################################

# Hint Project is set to 1 for projections

# Data vector
data <- list(Nyear=Nyear,Nage=Nage,M=M,Wght=Wght,SigCatch=SigCatch,SigCPUE=SigCPUE,Omega=Omega,
             Catch=Catch,CPUE=CPUE,Propn=Propn)
             
             
parameters <- list(dummy=0,LogN=LogN,Sel50=Sel50,Sel95=Sel95,LogFish=LogFish,logq=logq)

# When I was testing the code
#map<-list(LogN=rep(factor(NA),length(LogN)),Sel50=factor(NA),Sel95=factor(NA),
#                   LogFish=rep(factor(NA),length(LogFish)),logq=factor(NA))
# Estimate everything
map<-list(dummy=factor(NA))

#print(data)
#print(parameters)

################################################################################

model <- MakeADFun(data, parameters, DLL="LectC3",silent=T,map=map)

# test code - for checking for minimization
xx <- model$fn(model$env$last.par) #calling the function just to see what's happening
#print(model$report())
#cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")

#adding a bunch of parameters to make problem easier to solve
# Actual minimzation (with some "Bonus" parameters from nlminb)
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
best <- model$env$last.par.best
rep <- sdreport(model)
print(best)
print(rep)
print(model$report()$S)
print(model$report()$N)
print(model$report()$CPUEPred)
cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")
rep <- sdreport(model)
print(summary(rep))


