setwd("Project/Test/")

require(TMB)
compile("Test.cpp")
dyn.load(dynlib("Test"))

################################################################################

nyr <- scan("Test.dat",skip=1,n=1,quiet=T)
nage <- scan("Test.dat",skip=3,n=1,quiet=T)+1
M <- scan("Test.dat",skip=5,n=1,quiet=T)
data_cm_waa <- scan("Test.dat",skip=7,n=nage,quiet=T)
sig_catch <- scan("Test.dat",skip=9,n=1,quiet=T)
sig_cm_cpue <-scan("Test.dat",skip=11,n=1,quiet=T)  
omega_comp <-scan("Test.dat",skip=13,n=1,quiet=T)
data_catch <- matrix(scan("Test.dat",skip=15,n=nyr*3,quiet=T),ncol=3,byrow=T)[,2]
data_cm_cpue <- matrix(scan("Test.dat",skip=15,n=nyr*3,quiet=T),ncol=3,byrow=T)[,3]
data_cm_comp <- matrix(scan("Test.dat",skip=37,n=nyr*(nage+1),quiet=T),ncol=nage+1,byrow=T)[,-1]

################################################################################

logN <- scan("Test.pin",skip=2,n=nyr+nage-1,quiet=T)
cm_sel50 <- scan("Test.pin",skip=4,n=1,quiet=T)
cm_sel95 <- scan("Test.pin",skip=6,n=1,quiet=T)
logF <- scan("Test.pin",skip=8,n=nyr,quiet=T)
logq <- scan("Test.pin",skip=10,n=1,quiet=T)

################################################################################

# Data vector
data <- list(nyr=nyr,nage=nage,M=M,data_cm_waa=data_cm_waa,
             sig_catch=sig_catch,sig_cm_cpue=sig_cm_cpue,omega_comp=omega_comp,
             data_catch=data_catch,data_cm_cpue=data_cm_cpue,
             data_cm_comp=data_cm_comp)

parameters <- list(dummy=0,logN=logN,cm_sel50=cm_sel50,cm_sel95=cm_sel95,
                   logF=logF,logq=logq)

# When I was testing the code
#map<-list(logN=rep(factor(NA),length(logN)),cm_sel50=factor(NA),
    # cm_sel95=factor(NA),logF=rep(factor(NA),length(logF)),logq=factor(NA))
# Estimate everything
map<-list(dummy=factor(NA))

#print(data)
#print(parameters)

################################################################################

model <- MakeADFun(data, parameters, DLL="Test",silent=T,map=map)

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

# Model 1 ----

