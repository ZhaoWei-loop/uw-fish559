setwd("D:\\Research\\csiro\\Species18\\TMB Workshop\\QR\\")

#[a] treating recruitment as fixed or random 
#[b] allowing for time-varying q
#[c] include a PRI index of recruitment "properly"
#[d] estimating one or two sigmaCs (using the "map" function in TMB)

# Specify the case to run (0=fixed effects;1=random effects for rec;2=random effects for q and rec;3=random effects for rec+rec index)
TheCase <- 3

# Number of simulations  
Nsim <- 100                                             

# Number of years of data
Nyear <- 35                                            

# Set the variance parameters
SigC <- 60                                             # standard deviation of catch measurements
SigI <- 0.1                                            # CV of the index of abundance
SigQ <- 0.0                                            # variance of q
AR <- 0.9                                              # autocorrelation in q
SigmaR <- 0.5                                          # Recruitment variation

#==========================================================================================
#==========================================================================================

DoFit <- function(Case,Nyear,Nage,M,Wght,CatchW,CatchN,Effort,RecIndex)
{
 parameters <- list(dummy=0,F0=F0,LogRec1983=LogRec1983,LogRec=LogRec,logq=logq,LogSigCatch=LogSigCatch,
                    LogSigmaR=log(0.6),LogSigRecIndex=log(0.5),logqdev=logqdev,LogSigmaQ=log(0.1),LogitRho=0)
  
 # Default data and switches
 data <- list(Nyear=Nyear,Nage=Nage,M=M,Wght=Wght,CatchW=CatchW,CatchN=CatchN,RecIndex=RecIndex,
              Effort=Effort,RandomR=0,RandomQ=0,UseRecIndex=0)
  
 # fixed effects version
 if (Case==0)
  { 
   map<-list(dummy=factor(NA),LogSigmaR=factor(NA),LogSigRecIndex=factor(NA),
             logqdev=rep(factor(NA),Nyear-1),LogSigmaQ=factor(NA),LogitRho=factor(NA))
   model <- MakeADFun(data, parameters, random=NULL,DLL="LectH",silent=T,map=map)
  } 

 # Recruitment is a random effect
 if (Case==1)
  { 
   map<-list(dummy=factor(NA),LogSigRecIndex=factor(NA),
             logqdev=rep(factor(NA),Nyear-1),LogSigmaQ=factor(NA),LogitRho=factor(NA))
   data$RandomR = 1
   model <- MakeADFun(data, parameters, random="LogRec",DLL="LectH",silent=T,map=map)
  } 

 # Catchability is a random effect
 if (Case==2)
  { 
   map<-list(dummy=factor(NA),LogSigRecIndex=factor(NA))
   data$RandomR = 1
   data$RandomQ = 1
   model <- MakeADFun(data, parameters, random=c("LogRec","logqdev"),DLL="LectH",silent=T,map=map)
  } 
 
 # There is a recruitment index and recruitment is a random effect
 if (Case==3)
  { 
   map<-list(dummy=factor(NA),logqdev=rep(factor(NA),Nyear-1),LogSigmaQ=factor(NA),LogitRho=factor(NA))
   data$RandomR = 1
   data$UseRecIndex = 1
   model <- MakeADFun(data, parameters, random=c("LogRec"),DLL="LectH",silent=T,map=map)
  } 
  
 # Actual minimzation (with some "Bonus" parameters from nlminb)
 fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
 rep <- sdreport(model)
 #print(summary(rep))
 #AAA

 return(rep)
}

#==========================================================================================

require(TMB)
compile("LectH.cpp")
dyn.load(dynlib("LectH"))

################################################################################


Effort <- exp(rnorm(Nyear,0,0.5))                      # Effort

# set the seed
set.seed(999)

# basic biological parameters
Nage <- 6       # One less than the number of ages, i.e. for ages 1-20 set this to 19)
M <- 0.1
Wght <- c(1.0,1.0,1.0,2.0,2.0,2.0)

# Other parameters
TrueF0 <- 0.1
qbar <- 0.1
q <- rep(0,Nyear)
q[1] <- rnorm(1,0,SigQ)
for (Iyear in 2:Nyear) q[Iyear] <- AR*q[Iyear-1]+sqrt(1.0-AR^2)*rnorm(1,0,SigQ)

# Compute F, Z and exploitation rate by year`wa`
FF <- qbar*exp(q)*Effort
Z <- M + FF
Exploit <- FF/Z*(1.0-exp(-Z))

# Set the initial N
N <- matrix(0,ncol=Nage,nrow=Nyear+1)
R1983 <- 1000
N[1,1] <- R1983
for (Iage in 2:Nage) N[1,Iage] <- N[1,Iage-1]*exp(-M-TrueF0)
N[1,Nage] <- N[1,Nage]/(1.0-exp(-M-TrueF0))
N[2:(Nyear+1),1] <- R1983*exp(rnorm(Nyear,0,SigmaR))
for (Iyear in 1:Nyear)
 {
  for (Iage in 2:Nage) 
   N[Iyear+1,Iage] <- N[Iyear,Iage-1]*exp(-Z[Iyear]) 
  N[Iyear+1,Nage] <- N[Iyear+1,Nage] + N[Iyear,Nage]*exp(-Z[Iyear]) 
 }  

# store core outputs
Ntot <- apply(N,1,sum)[1:Nyear]
Bio <- apply(t(N)*Wght,2,sum)[1:Nyear]
CatchNT <- Ntot*Exploit
CatchWT <- Bio*Exploit
RecTrue <- N[1:Nyear,1]

################################################################################

F0 <- 0.2
LogRec1983<- log(2000)
logq <- log(0.2)
LogSigCatch <- log(5)
LogRec <- rep(LogRec1983,Nyear-1)
logqdev <- rep(0,Nyear-1)

################################################################################

# Data vector
data <- list(Nyear=Nyear,Nage=Nage,M=M,Wght=Wght,CatchW=CatchWT,CatchN=CatchNT,Effort=Effort,RandomR=0)
             
parameters <- list(dummy=0,F0=F0,LogRec1983=LogRec1983,LogRec=LogRec,logq=logq,LogSigCatch=LogSigCatch,
                   LogSigmaR=log(0.6),logqdev=logqdev,LogSigmaQ=log(0.1),LogitRho=0)

# When I was testing the code
map<-list(F0=factor(NA), LogRec1983=factor(NA),logq=factor(NA),LogSigCatch=factor(NA),LogSigmaR=factor(NA))

# Estimate everything
map<-list(dummy=factor(NA),LogSigmaR=factor(NA))

#print(data)
#print(parameters)
#print(map)

################################################################################

#model <- MakeADFun(data, parameters, DLL="LectH",silent=T,map=map)

# test code - for checking for minimization
#xx <- model$fn(model$env$last.par)

# Actual minimzation (with some "Bonus" parameters from nlminb)
#fit <- nlminb(model$par, model$fn, model$gr, random=NULL,control=list(eval.max=100000,iter.max=1000))
#best <- model$env$last.par.best
#rep <- sdreport(model)
#print(best)
#print(summary(rep))
#x <- summary(rep)

EstimatesB <- matrix(0,nrow=Nsim,ncol=Nyear)
EstimatesR <- matrix(0,nrow=Nsim,ncol=Nyear)
Estimatesq <- matrix(0,nrow=Nsim,ncol=Nyear)
Sigmas <- matrix(NA,nrow=Nsim,ncol=5)
seeds <- floor(runif(Nsim,100,10000))
for (Isim in 1:Nsim)
 {
  cat(Isim,seeds[Isim],"\n")
  set.seed(seeds[Isim])
  CatchN <- rnorm(Nyear,CatchNT,SigC) 
  CatchW <- rnorm(Nyear,CatchWT,SigC) 
  RecIndex <- exp(rnorm(Nyear,log(RecTrue),SigI))

  rep <- DoFit(TheCase,Nyear,Nage,M,Wght,CatchW,CatchN,Effort,RecIndex)

  x <- summary(rep)
  Ests <- x[rownames(x)=="Bio",1]
  EstimatesB[Isim,] <- Ests
  Ests <- x[rownames(x)=="Rec",1]
  EstimatesR[Isim,] <- Ests
  Ests <- x[rownames(x)=="qv",1]
  Estimatesq[Isim,] <- Ests
  Sigmas[Isim,1] <- x[rownames(x)=="SigmaR",1]
  Sigmas[Isim,2] <- x[rownames(x)=="SigCatch",1]
  Sigmas[Isim,3] <- x[rownames(x)=="SigmaQ",1]
  Sigmas[Isim,4] <- x[rownames(x)=="SigRecIndex",1]
}   
#print(Bio)
#print(EstimatesB)
#print(Estimatesq)
#print(EstimatesR)
#print(Sigmas)

stats <- rep(0,3)
for (Isim in 1:Nsim)
 {
  stats[1] <- stats[1] + sum(((EstimatesB[Isim,]-Bio)/Bio)^2) 
  stats[2] <- stats[2] + sum(((EstimatesR[Isim,]-RecTrue)/RecTrue)^2) 
  stats[3] <- stats[3] + sum(((Estimatesq[Isim,]-qbar*exp(q))/qbar*exp(q))^2) 
 }  
 stats[1] <- sqrt(stats[1]/(Nsim*Nyear))
 stats[2] <- sqrt(stats[2]/(Nsim*Nyear))
 stats[3] <- sqrt(stats[3]/(Nsim*Nyear))


par(mfrow=c(3,3))
ymax <- max(EstimatesB)*1.2
plot(1:Nyear,Bio,type="l",lwd=3,ylim=c(0,ymax),col="red")
for (Isim in 1:20) lines(1:Nyear,EstimatesB[Isim,],lty=2)  
lines(1:Nyear,Bio,col="red")

ymax <- max(EstimatesR)*1.2
plot(1:Nyear,RecTrue,type="l",lwd=3,ylim=c(0,ymax),col="red")
for (Isim in 1:20) lines(1:Nyear,EstimatesR[Isim,],lty=2)  
lines(1:Nyear,RecTrue,col="red")

ymax <- max(Estimatesq)*1.2
plot(1:Nyear,qbar*exp(q),type="l",lwd=3,ylim=c(0,ymax),col="red")
for (Isim in 1:20) lines(1:Nyear,Estimatesq[Isim,],lty=2)  
lines(1:Nyear,qbar*exp(q),col="red")

hist(Sigmas[,1],main="",xlab="SigmaR")
abline(v=SigmaR,lwd=3,col="red")
hist(Sigmas[,2],main="",xlab="Sigma Catch")
abline(v=SigC,lwd=3,col="red")
hist(Sigmas[,3],main="",xlab="Sigma q")
abline(v=SigQ,lwd=3,col="red")
hist(Sigmas[,4],main="",xlab="SigmaI")
abline(v=SigI,lwd=3,col="red")
cat(TheCase,round(stats*100,2))

