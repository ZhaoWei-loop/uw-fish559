
################################################################################

setwd("tmb-workshop-fhl/Day5")
require(TMB)

################################################################################

m <- scan("ex5.dat",skip=1,n=1,quiet=T)
TT <- scan("ex5.dat",skip=3,n=m,quiet=T)
Tmax <- scan("ex5.dat",skip=5,n=1,quiet=T)
B <- matrix(scan("ex5.dat",skip=7,n=m*Tmax,quiet=T),ncol=Tmax,byrow=T)
R <- matrix(scan("ex5.dat",skip=58,n=m*Tmax,quiet=T),ncol=Tmax,byrow=T)
Phi0 <- scan("ex5.dat",skip=109,n=m,quiet=T)

data <- list(m=m,TT=TT,Tmax=Tmax,B=B,R=R,Phi0=Phi0)

################################################################################
# Hint Nproj = 0 for basic estimation
# Basic estimation
################################################################################

# Set the parameters to initial values
mu<-0.2
log_tau<- 0.1
B0<-rep(500,m)
log_sigR<-rep(0.1,m)
eps<-rep(0.5,m)

parameters <- list(dummy=0,mu=mu,log_tau=log_tau,B0=B0,
                   log_sigR=log_sigR,eps=eps)

#print(data)
#print(parameters)

compile("Ex5Class.cpp")
dyn.load(dynlib("Ex5Class"))

# test code - for checking for minimization
map<-list(mu=factor(NA),log_tau=factor(NA),
          B0=rep(factor(NA),m),log_sigR=rep(factor(NA),m),
          eps=rep(factor(NA),m))
model <- MakeADFun(data, parameters, DLL="Ex5Class",silent=T,map=map)
model$report()

xx <- model$fn(model$env$last.par)
cat(xx,model$report()$obj_fun,"\n")
#AAA

# Now estimate everything
map<-list(dummy=factor(NA))
parameters
model <- MakeADFun(data, parameters, DLL="Ex5Class",silent=T,
                   map=map, random=c("eps"))

# Bounds on the parameters - needs to be equal to the number of ESTIMATED
# parameters not equal to the length of the total number of pars in the
# parameters vectors (dummy disappears) use:
length(model$par) # how many parameters for upper and lower bounds
parameters
lowbnd= c(-Inf, # mu
          -7, # tau
          rep(0.1, m), # B0
          rep(-7, m), # log sigR
          rep(-Inf, m)) # eps

uppbnd= c(Inf, # mu
          4, # tau
          rep(1000, m), # B0
          rep(4, m),
          rep(Inf, m)) # log sig

fit <- nlminb(model$par, model$fn, model$gr, 
              control=list(rel.tol=1e-12,
                           eval.max=100000,iter.max=10000),
              lower=lowbnd,upper=uppbnd)

best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(summary(rep))


