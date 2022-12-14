library(TMB)
setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Scripts")
parameters <- c(x=1)
m1 <- matrix(c(1,2,3,4),ncol=2,nrow=2)
m2 <- matrix(c(1,2,3,4,5,6),ncol=3,nrow=2)
print(m1 %*% m2)
data <- list(v1=c(1,2,2),v2=c(3,4,5),v3=c(5,6,7,8,9,10,11),m1=m1,m2=m2)
compile("LectB1.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("LectB1"))

################################################################################

model <- MakeADFun(data, parameters, DLL="LectB1",silent=T)
fit <- nlminb(model$par, model$fn, model$gr)

best <- model$env$last.par.best
rep <- sdreport(model)
print(model$report()$m4) #how to print different things you wanted to "report"
print(model$report()$m5)
print(model$report()$m6)
print(model$report()$m7)
print(t(m1))
print(m1)

print(best)
print(rep)