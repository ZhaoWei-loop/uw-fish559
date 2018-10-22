#devtools::install_github('colemonnahan/adnuts', build_vignettes=TRUE)
library(adnuts)
library(TMB)

# source("tools.R")
setwd("tmb-workshop-fhl/Day4")

# =================================================================================================================

MCMC <- function(model,covar,xxx,mult,Nreps=20000,Nburns=2000,Nevery=10)
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
  # Now consider mcmc sampling
  # post1 <- MCMC(model,VarCo,model$env$last.par.best,mult=1.0,Nreps=1000000,Nburns=100000,Nevery=1000)
  # post1a <- post1$Outs
  # save(post1a,file="post1.RData")

  # Now consider mcmc sampling (NUTS)
  library(adnuts)
  mcmcout <- sample_tmb(obj=model,
                        seeds=1:3,
                        init=list(list(model$env$last.par.best)), 
                        iter=10000,chains=1,algorithm="NUTS")
  
  post2 <- extract_samples(mcmcout)
  save(post2,file="post2.RData")
  # launch_shinytmb(mcmcout)
}  

#####----

# MCMCSumm <- function(file,best,Nyear,data,map,parameters)
# {
#   # Load the parameter file
#   load(file=file)
# 
#   # Set up for graphs
#   par(mfrow=c(2,2))
# 
#   Nsim <- length(post1a[,1])
#   Npar <- length(post1a[1,])
#  
#   # Good to check (compare with the sdreport)
#   for (II in 1:Npar)  
#    cat(II,best[II],mean(post1a[,II]),sd(post1a[,II]),"\n") 
# 
#   # Extract biomass and plot
#   Biomass <- matrix(0,nrow=Nsim,ncol=Nyear+1)
#   for (Isim in 1:Nsim)
#    { xx <- model$fn(post1a[Isim,]); Biomass[Isim,] <- model$report()$BioPred[1:(Nyear+1)]; }  
#   quant <- matrix(0,nrow=5,ncol=Nyear+1)
#   for (Iyear in 1:(Nyear+1))
#    quant[,Iyear] <- quantile(Biomass[,Iyear],probs=c(0.05,0.25,0.5,0.75,0.95))  
#   Years <- 1:(Nyear+1)
#   
#   # Plot of biomass posterior
#   ymax <- max(quant)
#   plot(Years,quant[3,],xlab="Year",ylab="Biomass",ylim=c(0,ymax),type='l')
#   xx <- c(Years,rev(Years))
#   yy <- c(quant[1,],rev(quant[5,]))
#   polygon(xx,yy,col="gray10")
#   xx <- c(Years,rev(Years))
#   yy <- c(quant[2,],rev(quant[4,]))
#   polygon(xx,yy,col="gray90")
#   lines(Years,quant[3,],lwd=3,lty=3)
#   
#   # Projection    
#   # ==========
#   
#   # hahaha - this is for you!
#   
# }  

##############################################################################
################################################################################

Nyear <- scan("ex4.dat",skip=1,n=1,quiet=T)
Nclass <- scan("ex4.dat",skip=3,n=1,quiet=T)
Length <- scan("ex4.dat",skip=5,n=Nclass,quiet=T)
Weight <- scan("ex4.dat",skip=7,n=Nclass,quiet=T)
X <- matrix(scan("ex4.dat",skip=9,n=Nclass*Nclass,quiet=T),ncol=Nclass,byrow=T)
M <- scan("ex4.dat",skip=19,n=1,quiet=T)
CWObs <- scan("ex4.dat",skip=22,n=Nyear,quiet=T)
CALObs <- matrix(scan("ex4.dat",skip=49,n=Nyear*Nclass,quiet=T),ncol=Nclass,byrow=T)
Neff <- scan("ex4.dat",skip=75,n=1,quiet=T)
BioIndex <- scan("ex4.dat",skip=78,n=Nyear,quiet=T)
BioSig <- scan("ex4.dat",skip=104,n=1,quiet=T)

# Selectivity parameter (convert to selectivity)
S50 <- scan("ex4.dat",skip=15,n=2,quiet=T)[1]
S95 <- scan("ex4.dat",skip=15,n=2,quiet=T)[2]
SS50 <- scan("ex4.dat",skip=17,n=2,quiet=T)[1]
SS95 <- scan("ex4.dat",skip=17,n=2,quiet=T)[2]
S <- rep(0,Nclass); SurveyS <- rep(0,Nclass)
for(Iclass in 1:Nclass)
 {
  S[Iclass] <- 1.0/(1+exp(-log(19.0)*(Length[Iclass]-S50)/(S95-S50)))
  SurveyS[Iclass] <- 1.0/(1+exp(-log(19.0)*(Length[Iclass]-SS50)/(SS95-SS50)))
 }  

# Normalize the catch-at-length data
for (Iyear in 1:Nyear)
 {
  Total <- 0
  for (Iclass in 1:Nclass) Total <- Total + CALObs[Iyear,Iclass]
  for (Iclass in 1:Nclass) CALObs[Iyear,Iclass] = CALObs[Iyear,Iclass]/Total
 }  

################################################################################
# Hint Nproj = 0 for basic estimation
# Basic estimation
################################################################################

# # Data
# Nyear <- readVec(file="EX4.DAT", "# Number of years")
# Nclass <- readVec(file="EX4.DAT", "# Number of length-classes")
# Length <- readVec(file="EX4.DAT", string="# Length-at-size")
# Weight <- readVec(file="EX4.DAT", string="# Weight-at-size")
# X <- readMat(file="EX4.DAT", string="# Size-transition matrix", Nclass)
# S <- readMat(file="EX4.DAT", string="# Size-transition matrix", Nclass)
# Set the parameters to initial values
LogRbar <- scan("ex4.pin",skip=1,n=1,quiet=T)
LogNinit <- scan("ex4.pin",skip=3,n=Nclass,quiet=T)
LogFullF <- scan("ex4.pin",skip=5,n=Nyear,quiet=T)
Eps <- scan("ex4.pin",skip=9,n=Nyear,quiet=T)
Fproj <- 0

# Data vector
data <- list(Nyear=Nyear,Nclass=Nclass,Length=Length,Weight=Weight,X=X,S=S,SurveyS=SurveyS,M=M,
             CWObs=CWObs,CALObs=CALObs,Neff=Neff,BioIndex=BioIndex,BioSig=BioSig,Nproj=0,Fproj=0,
             data=0,Jac_correction=0)
parameters <- list(dummy=0,LogRbar=LogRbar,LogNinit=LogNinit,LogFullF=LogFullF,Eps=Eps)
# When I was testing the code
map<-list(LogRbar=factor(NA),LogNinit=rep(factor(NA),Nclass),LogFullF=rep(factor(NA),Nyear),
          Eps=rep(factor(NA),Nyear))
# Estimate everything
map<-list(dummy=factor(NA))
map<-list()
#print(data)
#print(parameters)

compile("Ex4Class.cpp")
dyn.load(dynlib("Ex4Class"))

model <- MakeADFun(data, parameters, DLL="Ex4Class",silent=T,map=map,hessian=T)

# test code - for checking for minimization
xx <- model$fn(model$env$last.par)
#print(model$par)
#print(model$report())
#cat(model$report()$obj_fun,model$report()$LikeCatch,model$report()$LikeCAL,model$report()$LikeBio,model$report()$Penal,"\n")
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(rep)
cat(model$report()$obj_fun,model$report()$LikeCatch,model$report()$LikeCAL,model$report()$LikeBio,model$report()$Penal,"\n")
model$report()$rho

model$fn() # objective function


################################################################################
# Hint Nproj = 20 for projections
# Basic estimation
################################################################################

Nproj = 20
compile("Ex4Class.cpp")
dyn.load(dynlib("Ex4Class"))

# Set the parameters to initial values
LogRbar <- scan("ex4.pin",skip=1,n=1,quiet=T)
LogNinit <- scan("ex4.pin",skip=3,n=Nclass,quiet=T)
LogFullF <- scan("ex4.pin",skip=5,n=Nyear,quiet=T)
Eps <- scan("ex4.pin",skip=9,n=Nyear,quiet=T)
Eps <- c(Eps,rep(0,Nproj))

# Data vector
data <- list(Nyear=Nyear,Nclass=Nclass,Length=Length,Weight=Weight,X=X,S=S,SurveyS=SurveyS,M=M,
             CWObs=CWObs,CALObs=CALObs,Neff=Neff,BioIndex=BioIndex,BioSig=BioSig,Nproj=Nproj,Fproj=0,
             Jac_correction=1)
parameters <- list(dummy=0,LogRbar=LogRbar,LogNinit=LogNinit,LogFullF=LogFullF,Eps=Eps)
# When I was testing the code
map<-list(LogRbar=factor(NA),LogNinit=rep(factor(NA),Nclass),LogFullF=rep(factor(NA),Nyear),
          Eps=rep(factor(NA),Nyear))
# Estimate everything
map<-list(dummy=factor(NA))
map<-list()
#print(data)
#print(parameters)

# Note hessian=T so we get the Hessian matrix
model <- MakeADFun(data, parameters, DLL="Ex4Class",silent=T,map=map,hessian=T)

# test code - for checking for minimization (objective function)
xx <- model$fn(model$env$last.par)
#print(model$par)
#print(model$report())
#cat(model$report()$obj_fun,model$report()$LikeCatch,model$report()$LikeCAL,model$report()$LikeBio,model$report()$Penal,"\n")
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
best <- model$env$last.par.best
print(best)
cat(model$report()$obj_fun,model$report()$LikeCatch,model$report()$LikeCAL,model$report()$LikeBio,model$report()$Penal,"\n")
VarCo <- solve(model$he())
rep <- sdreport(model)
print(rep)
print(model$report()$obj_fun)
# Check for Hessian
print(sqrt(diag(VarCo)))
model$report()$nll
model$fn(best)#$best)
# Now consider mcmc sampling and provide a posterior for By
# ========================================================

#Switch correction to 1
data$Jac_correction <- 1
DoMCMCAll(model,VarCo)

# Now do the projections

# Load the parameter file
load(file="post2.RData")

# Set up for graphs
par(mfrow=c(2,2))

dim(post2)  

Nsim <- length(post2[,1])
Npar <- length(post2[1,])

# Rebuilding biomass (target)
target_biomass <- 1000
# Empty matrix to extract biomass from posteriors
Biomass <- matrix(0,nrow=Nsim,ncol=Nyear+Nproj)

# Vector of fishing probabilities
fs <- seq(0, 1.6, 0.015)

# Empty vector of probabilities
Probs <- rep(NA,length(fs))

#The aim of this step is to find the constant level of Fy, which if implemented
#over the next 20 years will result in By rebuilding to 1000t with 0.5
#probability by year 45, i.e. P(B45>1000t)=0.5.. For “bonus points” plot the
#future F as a function of P(B45>1000t).

for(i in 1:length(fs)) {

  data$Fproj <- fs[i]
  model <- MakeADFun(data, parameters, DLL="Ex4Class",
                     silent=T,map=map,hessian=T)
  for (j in 1:Nsim) {
    
    xx <- model$fn(post2[j,]) # Extract sample from posterior
    Biomass[j,] <- model$report()$BioPred
  }
  Probs[i] <- sum(Biomass[,Nyear+Nproj] > target_biomass) / Nsim 
}

plot(fs, Probs, type = "l", xlab = "Fishing mortality", ylab = "Probability < 1000 t")

# Bisection to find desire recovery probability
Fmin <- 0
Fmax <- 1

for (i in 1:20) {

  Fbar <- (Fmin + Fmax)/2.0
  data$Fproj <- Fbar
  model <- MakeADFun(data, parameters, DLL="Ex4Class",
                     silent=T,map=map,hessian=T)
  for (j in 1:Nsim) {
    xx <- model$fn(post2[j,]) # Extract sample from posterior
    Biomass[j,] <- model$report()$BioPred
  }
  Probs <- sum(Biomass[,Nyear+Nproj] > target_biomass) / Nsim
  if (Probs > 0.5) Fmin <- Fbar else Fmax <- Fbar
}
cat(Fbar,Probs,"\n")
points(Fbar,Probs,pch=16)
abline(h=0.5,lty=2)
abline(v=Fbar,lty=2)
Fbar

  
  # Good to check (compare with the sdreport)
  for (i in 1:Nsim)  
    cat(i,best[i],mean(post2[,i]),sd(post2[,i]),"\n") 
  
  # Extract biomass and plot
  Biomass <- matrix(0,nrow=Nsim,ncol=Nyear+Nproj)
  for (i in 1:Nsim) {
    xx <- model$fn(post2[i,]); 
    Biomass[i,] <- model$report()$BioPred[1:(Nyear+Nproj)]; 
  }  
  
  quant <- matrix(0,nrow=5,ncol=Nyear+Nproj)
  for (i in 1:(Nyear+Nproj)) {
    quant[,i] <- quantile(Biomass[,i],probs=c(0.05,0.25,0.5,0.75,0.95))  
    Years <- 1:(Nyear+Nproj)
  }
  
  # Plot of biomass posterior
  ymax <- max(quant)
  plot(Years,quant[3,],xlab="Year",ylab="Biomass",ylim=c(0,ymax),type='l')
  xx <- c(Years,rev(Years))
  yy <- c(quant[1,],rev(quant[5,]))
  polygon(xx,yy,col="gray10")
  xx <- c(Years,rev(Years))
  yy <- c(quant[2,],rev(quant[4,]))
  polygon(xx,yy,col="gray90")
  lines(Years,quant[3,],lwd=3,lty=3)
  
  # I was going back through the exercises and was wondering if you could
  # provide an example of using uniroot() in Ex. 4 from class to find the
  # P(B>1000)=0.5)