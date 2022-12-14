setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Examples/Ex1/Solution")


require(TMB)
#compile("Ex1.cpp")
compile("Ex1.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex1"))

################################################################################

Ndata <- scan("ex1.dat",skip=1,n=1)
print(Ndata)
dataD <- read.table("ex1.dat",skip=2,header=TRUE,nrow=Ndata)
Model_type = 2 #starting w model type 2 i guess
data <- list(Age=dataD$Age,Length=dataD$Length,Ndata=Ndata,Model_type=Model_type) #supplying data
if (Model_type==1) map <- list(LogKappa=factor(NA),t0=factor(NA)) #cutting out params that don't apply to those eq
if (Model_type==2) map <- list(Loga50=factor(NA),LogDelta=factor(NA))
parameters <- list(LogLinf=4.78,Loga50=2.3,LogDelta=2.1,LogKappa=-3,t0=0, LogSigma=0) #initializing all

################################################################################

#here the models are being done individually, no if statements on the map arguments
Model_type = 1
data <- list(Age=dataD$Age,Length=dataD$Length,Ndata=Ndata,Model_type=Model_type)
map <- list(LogKappa=factor(NA),t0=factor(NA))
parameters <- list(LogLinf=4.78,Loga50=2.3,LogDelta=2.1,LogKappa=-3,t0=0, LogSigma=0)
model <- MakeADFun(data, parameters, DLL="Ex1",silent=T,map=map)
fit <- nlminb(model$par, model$fn, model$gr)
print(fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep)) #summarizes the whole model w standard errors
ThePred1 <- model$report()$PredY
ThePred <- model$report()$Pred #trying something

Model_type = 2
data <- list(Age=dataD$Age,Length=dataD$Length,Ndata=Ndata,Model_type=Model_type)
map <- list(Loga50=factor(NA),LogDelta=factor(NA))
parameters <- list(LogLinf=4.78,Loga50=2.3,LogDelta=2.1,LogKappa=-3,t0=0, LogSigma=0)
model <- MakeADFun(data, parameters, DLL="Ex1",silent=T,map=map)
fit <- nlminb(model$par, model$fn, model$gr)
print(fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))
ThePred2 <- model$report()$PredY
print(ThePred2)
print(length(ThePred2))

par(mfrow=c(2,2))
plot(data$Age,data$Length,xlab="Age",ylab="Length",pch=16)
lines(1:20,ThePred1,lty=1)
lines(1:20,ThePred2,lty=2)
legend("topleft",legend=c("Logistic","Von Bertlanffy"),lty=1:2)

print(ThePred2)