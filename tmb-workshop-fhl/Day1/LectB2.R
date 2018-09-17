setwd("D:\\courses\\FISH 559_18\\TMB Workshop\\Lecture Examples\\")
data <- read.table("LectB2.dat", header=TRUE)
parameters <- list(b0=0, b1=0, logSigma=0)

# Read in data
# Specify parameters and starting values
# Create a dll
# 


require(TMB)
compile("LectB2.cpp")
dyn.load(dynlib("LectB2"))

################################################################################

# DLL = tells you which model to use
# silent = TRUE prevents you from getting a lot of unwanted output
model <- MakeADFun(data, parameters, DLL="LectB2",silent=T) 
fit <- nlminb(model$par, model$fn, model$gr)

best <- model$env$last.par.best
rep <- sdreport(model)

print(best)
print(rep)
