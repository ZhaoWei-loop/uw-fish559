setwd("tmb-workshop-fhl/Day2")

require(TMB)

compile("ex2v2.cpp")
dyn.load(dynlib("ex2v2"))

################################################################################

nprey <- scan("EX2.dat",skip=1,n=1)
Ndata <- scan("EX2.dat",skip=3,n=1)
data <- read.table("EX2.dat",skip=4,header=TRUE,nrow=Ndata)
colnames(data) <- c("Predator", "Prey1", "Prey2",
                    "Prey3", "C1", "C2", "C3")

Model_type = 4

mod_dat <- list(Predator = data[,1],
                Prey = as.matrix(data[,2:4]),
                Consump = as.matrix(data[,5:7]))

# Add Model_type to our list of data
mod_dat$Model_type <- Model_type
mod_dat$Ndata <- Ndata

# Parameter list
parameters <- list(logalpha = rep(0, 3), logbeta = rep(0, 3), gamma = rep(2.1, 3),
                   loglambda = rep(0, 3))

# Create maps for different model types


if (Model_type==1) map <- list(logbeta = factor(rep(NA, 3)), gamma = factor(rep(NA, 3)), 
                               loglambda = factor(rep(NA, 3)))
if (Model_type==2) map <- list(gamma = factor(rep(NA, 3)), loglambda = factor(rep(NA, 3)))
if (Model_type==3) map <- list(loglambda = factor(rep(NA, 3)))
if (Model_type==4) map <- list(gamma = factor(rep(NA, 3)))
################################################################################

model <- MakeADFun(mod_dat, parameters, DLL="ex2v2",map=map)
fit <- nlminb(model$par, model$fn, model$gr)
best <- model$env$last.par.best
rep <- sdreport(model)
print(rep)

# Model 4 predictions
summary(rep)
Pred <- model$report()

#Predict consumption from 0 to 3.5
length(Pred$Pred)
plot(data$C1, data$Predator)
lines(data$C1, Pred$Pred)

plot(data$C2, data$Predator)
plot(data$C3, data$Predator)
