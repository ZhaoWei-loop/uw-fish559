setwd("D:\\courses\\FISH 559_18\\TMB Workshop\\Lecture Examples\\")




#==========================================================================================

require(TMB)
compile("LectD.cpp")
dyn.load(dynlib("LectD"))

################################################################################

Nyear <- scan("LectD.dat",skip=1,n=1,quiet=T)
Nage <- scan("LectD.dat",skip=3,n=1,quiet=T)+1
M <- scan("LectD.dat",skip=5,n=1,quiet=T)
Wght <- scan("LectD.dat",skip=7,n=Nage,quiet=T)
SigCatch <- scan("LectD.dat",skip=9,n=1,quiet=T)
SigCPUE <-scan("LectD.dat",skip=11,n=1,quiet=T)  
Omega <-scan("LectD.dat",skip=13,n=1,quiet=T)
Catch <- matrix(scan("LectD.dat",skip=15,n=Nyear*3,quiet=T),ncol=3,byrow=T)[,2]
CPUE <- matrix(scan("LectD.dat",skip=15,n=Nyear*3,quiet=T),ncol=3,byrow=T)[,3]
Propn <- matrix(scan("LectD.dat",skip=37,n=Nyear*(Nage+1),quiet=T),ncol=Nage+1,byrow=T)[,-1]

################################################################################

LogN <- scan("LectD.pin",skip=2,n=Nyear+Nage-1,quiet=T)
Sel50 <- scan("LectD.pin",skip=4,n=1,quiet=T)
Sel95 <- scan("LectD.pin",skip=6,n=1,quiet=T)
LogFish <- scan("LectD.pin",skip=8,n=Nyear,quiet=T)
logq <- scan("LectD.pin",skip=10,n=1,quiet=T)

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

model <- MakeADFun(data, parameters, DLL="LectD",silent=T,map=map)

# test code - for checking for minimization
xx <- model$fn(model$env$last.par)
#print(model$report())
#cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")

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


