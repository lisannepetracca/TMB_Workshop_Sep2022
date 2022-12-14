#devtools::install_github('colemonnahan/adnuts', build_vignettes=TRUE)

# =================================================================================================================
# =================================================================================================================

Leapfrog <- function(Theta,rr,eps)
 {
  rtilda <- rr - eps/2*as.numeric(model$gr(Theta))
  ThetaTilda <- Theta + eps*rtilda
  rtilda2 <- rtilda - eps/2*as.numeric(model$gr(ThetaTilda)) 
  Outs <- NULL
  Outs$ThetaTilda <- ThetaTilda
  Outs$rtilda <- rtilda2
  return(Outs)
 }  

# =================================================================================================================
# =================================================================================================================

BuildTree <- function(model,Theta,ru,ThetaMinus,ThetaPlus,rminus,rplus,u,vj,j,eps)
{
 DeltaMax <- 1000
 logu <- log(u) 
 if (j==0)
  {
   cat("j=0",Theta,"\n")
   Outs <- Leapfrog(Theta,ru,eps)
   rtilda <- Outs$rtilda; ThetaTilda <- Outs$ThetaTilda
   LogL <- -1*model$fn(ThetaTilda)
   if (logu < LogL - sum(rtilda*rtilda)/2.0)
    Cd <- c(ThetaTilda,rtilda) 
   else
    Cd <- NULL
   Sdash <- ifelse(logu < DeltaMax + LogL - sum(rtilda*rtilda)/2,1,0)
   rminus <- rtilda
   ThetaMinus <- ThetaTilda
   rplus <- rtilda
   ThetaPlus <- ThetaTilda
  } 
 else
  {
   Outs <- BuildTree(model,Theta,ru,ThetaMinus,ThetaPlus,rminus,rplus,u,vj,j-1,eps)
   ThetaMinus <- Outs$ThetaMinus; rminus <- Outs$rminus; ThetaPlus <- Outs$ThetaPlus; rplus <- Outs$rplus
   Cd <- Outs$Cd; Sdash <- Outs$Sdash
   if (vj < 0)
    { 
     cat("vj negative",Theta,"\n")
     Outs <- BuildTree(model,ThetaMinus,rminus,ThetaMinus,ThetaPlus,rminus,rplus,u,vj,j-1,eps)
     ThetaMinus <- Outs$ThetaMinus; rminus <- Outs$rminus; Cdd <- Outs$Cd; Sdash2 <- Outs$Sdash
     
    } 
   else
    {
     cat("vj postive",Theta,"\n")
     Outs <- BuildTree(model,ThetaPlus,rplus,ThetaMinus,ThetaPlus,rminus,rplus,u,vj,j-1,eps)
     ThetaPlus <- Outs$ThetaPlus; rplus <- Outs$rplus; Cdd <- Outs$Cd; Sdash2 <- Outs$Sdash
    }
   Term1 <- sum((ThetaMinus-ThetaPlus)*rminus)
   Term2 <- sum((ThetaMinus-ThetaPlus)*rplus)
   Sdash <- Sdash*Sdash2*ifelse(Term1>=0,1,0)*ifelse(Term2>=0,1,0)
   Cd <- rbind(Cd,Cdd);
  } 

 Out2 <- NULL
 Out2$ThetaMinus <- ThetaMinus
 Out2$rminus <- rminus
 Out2$ThetaPlus <- ThetaPlus
 Out2$rplus <- rplus
 Out2$Cd <- Cd
 Out2$Sdash <- Sdash
 return(Out2)
}  

# ------------------------------------------------------------------------------------------

MCMCNUTS <- function(model,xxx,Nreps=22000,Nburns=2000,Nevery=10,eps=0.01)
{
  print("Doing Simple NUTS")  
  Npars <- length(xxx)
  Nout <- Nreps
  set.seed(666)
  
  # Initialize
  Theta <- as.numeric(xxx)
  gold <- exp(-1*model$fn(Theta))
  Outputs <- matrix(0,nrow=Nout,ncol=Npars)
  ObjFn <- rep(NA,length=Nout)
  
  Iout <- 0
  Itst <- 0
  for (Irep in 1:Nreps)
   { 
    C <- matrix(0,nrow=1,ncol=(Npars+Npars))
    r0 <- rnorm(Npars,0,1)
    uu <- runif(1,0,gold*exp(-sum(r0*r0)/2))
    cat("rep = ",Irep,"\n")

    # Initialize 
    ThetaMinus <- Theta; ThetaPlus <- Theta; rminus <- r0; rplus <- r0; j <- 0; C[1,] <- c(Theta,r0); s <- 1
    while (s==1)
     {
      vj <- runif(1,-1,1)
      if (vj < 0)
       {
        StepNeg <- BuildTree(model,ThetaMinus,rminus,ThetaMinus,ThetaPlus,rminus,rplus,uu,vj,j,eps)        
        ThetaMinus <- StepNeg$ThetaMinus; rminus <- StepNeg$rminus; Cd <- StepNeg$Cd; Sdash <- StepNeg$Sdash
      }
      else
       {
        StepPos <- BuildTree(model,ThetaPlus,rplus,ThetaMinus,ThetaPlus,rminus,rplus,uu,vj,j,eps)        
        ThetaPlus <- StepPos$ThetaPlus; rplus <- StepPos$rplus; Cd <- StepPos$Cd; Sdash <- StepPos$Sdash
       }  
      if (Sdash == 1) C <- rbind(C,Cd)
      Term1 <- sum((ThetaMinus-ThetaPlus)*rminus)
      Term2 <- sum((ThetaMinus-ThetaPlus)*rplus)
      #cat(Sdash,Term1,Term2,"\n")
      s <- Sdash*ifelse(Term1>=0,1,0)*ifelse(Term2>=0,1,0)
      j <- j + 1
     }  # End s loop 
    #print(C)
    Index <- floor(runif(1,0,length(C[,1]))+1)
    Theta <- C[Index,1:Npars]
    gold <- exp(-1*model$fn(Theta))
    if (length(C[,1])>1000)
     {  
      Itst <- Itst + 1
      if (Itst==1) par(mfrow=c(3,3))
      Index <- length(C[,1])
      plot(C[,1],C[,2],type='b',pch=2)
      points(C[1,1],C[1,2],pch=16)
      points(C[Index,1],C[Index,2],pch=18,cex=4)
      #if (Itst==9) AAA
     }  

    # end of current replicate
    if (Irep %% Nevery == 0)
      if (Irep > Nburns)
      { 
        Iout <- Iout + 1 
        Outputs[Iout,] <- Theta 
        ObjFn[Iout] <- log(gold)
      }
    if (Irep %% 1000 == 0) print(Irep)  
  } 
  
  Outs <- list()
  Outs$Outs <- Outputs
  Outs$Outs <- Outs$Outs[1:Iout,]
  Outs$ObjFn <- ObjFn[1:Iout]
  return(Outs)
}  
  
# =================================================================================================================

MCMCHMC <- function(model,xxx,Nreps=22000,Nburns=2000,Nevery=10,eps=0.005,L=5)
{
  Npars <- length(xxx)
  Nout <- Nreps
  set.seed(666)

  # Counters
  ichk <- 0
  jchk <- 0

    # Initialize
  Theta <- as.numeric(xxx)
  gold <- exp(-1*model$fn(Theta))
  Outputs <- matrix(0,nrow=Nout,ncol=Npars)
  ObjFn <- rep(NA,length=Nout)
  
  Iout <- 0
  for (Irep in 1:Nreps)
   { 
    r0 <- rnorm(Npars,0,1)
    ThetaTilda <- Theta
    ThetaOld <- Theta
    rtilda <- r0
    for (II in 1:L)
     {
      Outs <- Leapfrog(ThetaTilda,rtilda,eps)
      rtilda <- Outs$rtilda; ThetaTilda <- Outs$ThetaTilda
     }  
    gold1 <- exp(-1*model$fn(ThetaTilda))
    Ratio <- gold1*exp(-sum(rtilda*rtilda)/2.0)/(gold*exp(-sum(r0*r0)/2.0))
    RND <-runif(1,0,1)
    if (Ratio > RND)
    { gold <- gold1; Theta <- ThetaTilda;}
    else
    { Theta <- ThetaOld; } 

    # end of current replicate
    if (Irep %% Nevery == 0)
      if (Irep > Nburns)
      { 
        Iout <- Iout + 1 
        Outputs[Iout,] <- Theta 
        ObjFn[Iout] <- log(gold)
      }
    if (Irep %% 1000 == 0) print(Irep)  
   } 
  
  Outs <- list()
  Outs$Outs <- Outputs
  Outs$Outs <- Outs$Outs[1:Iout,]
  Outs$ObjFn <- ObjFn[1:Iout]
  return(Outs)
}  

# =================================================================================================================

MCMCNaive <- function(model,covar,xxx,mult,Nreps=22000,Nburns=2000,Nevery=10)
{
  library(mvtnorm)
  
  Npars <- length(xxx)
  Nout <- Nreps
  Outputs <- matrix(0,nrow=Nout,ncol=Npars)
  ObjFn <- rep(NA,length=Nout)
  
  # Counters
  ichk <- 0
  jchk <- 0
  
  # Initialize
  Theta <- xxx
  DeltaInit <- 0.1
  Delta <- DeltaInit * xxx
  print(Theta)
  gold <- exp(-1*model$fn(Theta))
  print("gold")
  print(log(gold))
  
  Iout <- 0
  for (Irep in 1:Nreps)
  {
    
    # do for each parameter
    asave <- Theta
    Theta <- rmvnorm(1, mean=asave, sigma=covar*mult)
    g1 <- exp(-1*model$fn(Theta))
    #print(log(g1))
    
    RND <-runif(1,0,1)
    if (g1/gold > RND)
    { gold <- g1; ichk <- ichk+1 }
    else
    { Theta <- asave; jchk <- jchk + 1 } 
    
    if (Irep <= Nburns) 
      if (Irep %% Nevery == 0)
      {
        if (ichk > jchk)
          mult <- mult * 1.05
        else 
          mult <- mult * 0.95
        ichk <- 0
        jchk <- 0
      }  
    
    # end of current replicate
    if (Irep %% Nevery == 0)
      if (Irep > Nburns)
      { 
        Iout <- Iout + 1 
        Outputs[Iout,] <- Theta 
        ObjFn[Iout] <- log(gold)
      }
    if (Irep %% 1000 == 0) print(Irep)  
    
  }
  
  Outs <- list()
  Outs$Outs <- Outputs
  Outs$Outs <- Outs$Outs[1:Iout,]
  Outs$ObjFn <- ObjFn[1:Iout]
  return(Outs)
}

################################################################################

DoMCMCAll <- function(model,VarCo)
 {  
  windows(height=10,width=10)
  par(mfrow=c(4,4),mar=c(4,3,1,1))
  
  
  # Now consider NUts sampling
  post0 <- MCMCNUTS(model,model$env$last.par.best,Nreps=20000,Nburns=0,Nevery=20)
  post1a <- post0$Outs
  print(head(post1a))
  Nrow <- length(post1a[,1])  
  xx <- seq(1,Nrow)
  plot(xx,exp(post1a[,1]),xlab="Cycle number",ylab="Parameter 1",pch=16,ylim=c(35,42))
  plot(xx,exp(post1a[,2]),xlab="Cycle number",ylab="Parameter 2",pch=16,ylim=c(65,74))
  plot(xx,exp(post1a[,3]),xlab="Cycle number",ylab="Parameter 3",pch=16,ylim=c(0.1,0.4))
  plot(xx,-1*post0$ObjFn,xlab="Cycle number",ylab="Post Density",pch=16)
  dataout<-matrix(c(exp(post1a[,1]),exp(post1a[,2]),exp(post1a[,3]),-1*post0$ObjFn),nrow=Nrow)
  pairs(dataout,labels=c("A50","A95","Sigma","Post Density"),pch=16)
  
  
  # Now consider hmc sampling
  post0 <- MCMCHMC(model,model$env$last.par.best,Nreps=20000,Nburns=1,Nevery=20)
  post1a <- post0$Outs
  print(head(post1a))
  Nrow <- length(post1a[,1])  
  xx <- seq(1,Nrow)
  plot(xx,exp(post1a[,1]),xlab="Cycle number",ylab="Parameter 1",pch=16,ylim=c(35,42))
  plot(xx,exp(post1a[,2]),xlab="Cycle number",ylab="Parameter 2",pch=16,ylim=c(65,74))
  plot(xx,exp(post1a[,3]),xlab="Cycle number",ylab="Parameter 3",pch=16,ylim=c(0.1,0.4))
  plot(xx,-1*post0$ObjFn,xlab="Cycle number",ylab="Post Density",pch=16)
  dataout<-matrix(c(exp(post1a[,1]),exp(post1a[,2]),exp(post1a[,3]),-1*post0$ObjFn),nrow=Nrow)
  pairs(dataout,labels=c("A50","A95","Sigma","Post Density"),pch=16)
  AAA
  
  
  # Now consider mcmc sampling
  post1 <- MCMCNaive(model,VarCo,model$env$last.par.best,mult=1.0,Nreps=20000,Nburns=1,Nevery=20)
  post1a <- post1$Outs
  print(head(post1a))
  Nrow <- length(post1a[,1])  
  xx <- seq(1,Nrow)
  plot(xx,exp(post1a[,1]),xlab="Cycle number",ylab="Parameter 1",pch=16,ylim=c(35,42))
  plot(xx,exp(post1a[,2]),xlab="Cycle number",ylab="Parameter 2",pch=16,ylim=c(65,74))
  plot(xx,exp(post1a[,3]),xlab="Cycle number",ylab="Parameter 3",pch=16,ylim=c(0.1,0.4))
  plot(xx,-1*post1$ObjFn,xlab="Cycle number",ylab="Post Density",pch=16)
  dataout<-matrix(c(exp(post1a[,1]),exp(post1a[,2]),exp(post1a[,3]),-1*post1$ObjFn),nrow=Nrow)
  #pairs(dataout,labels=c("A50","A95","Sigma","Post Density"),pch=16)
  

  # Now consider mcmc sampling (NUTS)
  library(adnuts)
  mcmcout <- sample_tmb(obj=model,seeds=1:3,init=list(list(model$env$last.par.best)), iter=2000,chains=1,algorithm="NUTS")
  post2 <- extract_samples(mcmcout)
  print(head(post2))
  Nrow <- length(post2[,1])  
  xx <- seq(1,Nrow)
  plot(xx,exp(post2[,1]),xlab="Cycle number",ylab="Parameter 1",pch=16,ylim=c(35,42))
  plot(xx,exp(post2[,2]),xlab="Cycle number",ylab="Parameter 2",pch=16,ylim=c(65,74))
  plot(xx,exp(post2[,3]),xlab="Cycle number",ylab="Parameter 3",pch=16,ylim=c(0.1,0.4))
  dataout<-matrix(c(exp(post2[,1]),exp(post2[,2]),exp(post2[,3])),nrow=Nrow)
  #pairs(dataout,labels=c("A50","A95","Sigma"),pch=16)
}  

################################################################################

setwd("C:\\courses\\FISH 559_22\\TMB Workshop\\Lecture Examples\\")

require(TMB)

################################################################################

TheData<-scan(file="lectF.txt",what=list(Length=0,Prob=0))
data <- list(Length=TheData$Length,Prob=TheData$Prob,Prob2=log(TheData$Prob/(1-TheData$Prob)))
parameters = list(LogA50=log(40),LogA95=log(80),LogSigma=log(5))

print(data)
print(parameters)

compile("LectF2.cpp")
dyn.load(dynlib("LectF2"))

map <- list(LogSigma=factor(NA))
model <- MakeADFun(data, parameters, DLL="LectF2",silent=T,hessian=T,map=map)
map <- list(LogSigma=factor(NA))
model <- MakeADFun(data, parameters, DLL="LectF2",silent=T,hessian=T)
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
best <- model$env$last.par.best
rep <- sdreport(model)
print(rep)
VarCo <- solve(model$he())

# Now consider NUts sampling
post0 <- MCMCNUTS(model,model$env$last.par.best,Nreps=20000,Nburns=0,Nevery=20)


# Now consider mcmc sampling
# ==========================
DoMCMCAll(model,VarCo)


  