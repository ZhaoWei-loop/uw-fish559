source("tools.R")

setwd("tmb-workshop-fhl/Day2")

N <- readVec(file = "ex1.dat", "Number of data points")
dat <- readMat(file = "ex1.dat", "Age    Length", nrow = N)

dat <- data.frame(dat)
colnames(dat) <- c("age", "len")

library(dplyr)
dat %>% 
  filter(len > 0) -> dat

plot(dat$age ~ dat$len)

data <- list(model = 1, len = dat$len, age = dat$age)
parameters <- list(logLinf = 5, logSigma = 0, logK = 0, loga0 = 1, loga50 = 2, logdelta = 3)

# map <- list(logK = factor(NA), loga0 = factor(NA), loga50 = factor(NA), logdelta = factor(NA))
map <- list(loga50 = factor(NA), logdelta = factor(NA))
# map <- list(logK = factor(NA), loga0 = factor(NA), loga50 = factor(NA), logdelta = factor(NA))
print(parameters)


require(TMB)
compile("ex1.cpp")
dyn.load(dynlib("ex1"))

model <- MakeADFun(data, parameters, DLL="ex1",
                   map = map)#,
                   # increase time to minimize and increase tolerance (get closer to minimum)
                   #control=list(eval.max=10000,iter.max=1000,rel.tol=1e-15),silent=T)
print(attributes(model))

fit <- nlminb(model$par, model$fn, model$gr)
for (i in 1:3)
  fit <- nlminb(model$env$last.par.best, model$fn, model$gr)

rep <- sdreport(model)

summary(rep)
# Sumamrize ALL
print(summary(rep,p.value=T))
# Restrict what comes out to the fixed parameters only
print(summary(rep,select="fixed",p.value=F))
