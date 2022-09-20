setwd("G:/My Drive/GitHub/TMB_Workshop_Sep2022/Examples/Ex2")


require(TMB)
compile("Ex2.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex2")) #it compiles but then crashes at the model step down below

################################################################################


dataD <- read.table("ex2.dat",skip=4,header=TRUE)
dim(dataD)
#Model_type = 1
Ndata <- nrow(dataD) 
Nprey <- 3
data <- list(Predator=dataD$Predator,Consump=as.matrix(dataD[,5:7]),Ndata=Ndata,Nprey=Nprey)#,Model_type=Model_type)
#if (Model_type==1) map <- list(LogKappa=factor(NA),t0=factor(NA))
#if (Model_type==2) map <- list(Loga50=factor(NA),LogDelta=factor(NA))
parameters <- list(alpha=c(1,2,3))
model <- MakeADFun(data, parameters, DLL="Ex2",silent=T)#,map=map) #ok so it crashes here and that is sad
fit <- nlminb(model$par, model$fn, model$gr)
################################################################################


