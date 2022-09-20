

setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Examples/Ex3")

dataD <- read.table("ex3.dat",skip=3,header=TRUE)
dataD
dim(dataD)
colnames(dataD) <- c("year", "Google", "WOS")
head(dataD)

require(TMB)
#compile("Ex1.cpp")
compile("Ex3.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex3"))

################################################################################

Ndata <- nrow(dataD)
data <- list(Google=dataD$Google, Ndata=Ndata)#,Model_type=Model_type)

parameters <- list(logLambda=5) #only one parameter
model <- MakeADFun(data, parameters, DLL="Ex3",silent=T)#,map=map)
fit <- nlminb(model$par, model$fn, model$gr)

par(mfrow=c(1,1))
plot(dataD$year,dataD$Google,xlab="Year",ylab="Citations",pch=16) #mapping year vs citations
ThePred1 <- model$report()$lambda #retrieving predictions from the model
lines(1982:2014,rep(ThePred1,Ndata),lty=2) #plotting them w data

NLL <- fit$objective #this returns NLL
(AICc <- (2*(NLL) + (2*1)) + ((2*1*1)+(2*1))/(33-1-1)) #this calculates AICc by hand; k is 1 (lambda), N is 33

################################################################################
