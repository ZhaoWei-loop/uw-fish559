setwd("tmb-workshop-fhl/Day2")

require(TMB)


################################################################################

nprey <- scan("EX2.dat",skip=1,n=1)
Ndata <- scan("EX2.dat",skip=3,n=1)
data <- read.table("EX2.dat",skip=4,header=TRUE,nrow=Ndata)
colnames(data) <- c("Predator", "Prey1", "Prey2",
                    "Prey3", "C1", "C2", "C3")

# INPUTS

Species = 1
Model_type = 4

if (Species==1) mod_dat <- list(Predator = data$Predator,
                                Prey = data$Prey1,
                                Consump = data$C1)

if (Species==2) mod_dat <- list(Predator = data$Predator,
                                Prey = data$Prey2,
                                Consump = data$C2)

if (Species==3) mod_dat <- list(Predator = data$Predator,
                                Prey = data$Prey3,
                                Consump = data$C3)
# Add Model_type to our list of data
mod_dat$Model_type <- Model_type
mod_dat$Ndata <- Ndata

# Create maps for different model types
if (Model_type==1) map <- list(logbeta = factor(NA), gamma = factor(NA), loglambda = factor(NA))
if (Model_type==2) map <- list(gamma = factor(NA), loglambda = factor(NA))
if (Model_type==3) map <- list(loglambda = factor(NA))
if (Model_type==4) map <- list(gamma = factor(NA))

# Parameter list
parameters <- list(logalpha = 0, logbeta = 2, gamma = 2.1, loglambda = 1)

################################################################################


compile("ex2.cpp")
dyn.load(dynlib("ex2"))
model <- MakeADFun(mod_dat, parameters, DLL="ex2",
                   silent=TRUE,map=map)
fit <- nlminb(model$par, model$fn, model$gr)
best <- model$env$last.par.best
rep <- sdreport(model)
print(rep)

# Model 4 predictions
summary(rep)
model$report()

max(data$C1)
max(data$C2)
max(data$C3)
#Predict consumption from 0 to 3.5

plot(data$C1, data$Predator)
plot(data$C2, data$Predator)
plot(data$C3, data$Predator)

Model_type = 2
data <- list(Age=dataD$Age,Length=dataD$Length,Ndata=Ndata,Model_type=Model_type)
map <- list(Loga50=factor(NA),LogDelta=factor(NA))
parameters <- list(LogLinf=4.78,Loga50=2.3,LogDelta=2.1,LogKappa=-3,t0=0, LogSigma=0)
model <- MakeADFun(data, parameters, DLL="Ex1",silent=T,map=map)
fit <- nlminb(model$par, model$fn, model$gr)
best <- model$env$last.par.best
rep <- sdreport(model)
print(rep)
ThePred2 <- model$report()$PredY
print(ThePred2)
print(length(ThePred2))

par(mfrow=c(2,2))
plot(data$Age,data$Length,xlab="Age",ylab="Length",pch=16)
lines(1:20,ThePred1,lty=1)
lines(1:20,ThePred2,lty=2)
legend("topleft",legend=c("Logistic","Von Bertlanffy"),lty=1:2)

print(ThePred2)