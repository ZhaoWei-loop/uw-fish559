# setwd("D:\\courses\\FISH 559_18\\TMB Workshop\\In Class Examples\\")
setwd("~/uw-fish559/tmb-workshop-fhl/Day1")
data <- list(x=rivers)
# initial values. using log sigma, which you never want ot be negative
parameters <- list(mu=0,logSigma=0) 

require(TMB)
compile('LectB1.cpp')
#load dyamic link library
dyn.load(dynlib('LectB1'))

##################
# creates the functions that are used for minimization
model <- MakeADFun(data,parameters,silent=T)
# minimization
fit   <- nlminb(model$par, model$fn, model$gr)
# gradient is the derivative matrix
rep   <- sdreport(model)
print(summary(rep))
