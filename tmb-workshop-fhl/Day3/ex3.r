setwd("tmb-workshop-fhl/Day3")
require(TMB)

# Compile model
compile("ex3.cpp")
dyn.load(dynlib("ex3"))

# read in data
nyears <- scan("ex3.dat",skip=1,n=1)
njourn <- 1
data <- read.table("ex3.dat",skip=5,nrow=nyears)
colnames(data) <- c("year", "google", "web")
data$year <- data$year - min(data$year) + 1

# plot data
par(mfrow=c(1, 2))
plot(data$year, data$google, type = 'p')
plot(data$year, data$web, type = 'p')

#put data in lists
mod_data <- list(model = 1, nyears = nyears, njourn = njourn, 
                 data = as.matrix(data)[,1:2])

#parameter section
parameters <- list(b0 = rep(0, njourn), log_b1 = rep(1, njourn), 
                   amp = rep(0.5, njourn), phase = rep(0, njourn), 
                   log_period = rep(log(4), njourn))

# map <- list(log_period = factor(rep(NA, njourn)))

model <- MakeADFun(mod_data, parameters, DLL="ex3",
                   # map = map,
                   silent = T, hessian = T)

fit <- nlminb(model$par, model$fn, model$gr,
              control=list(eval.max=100000,iter.max=1000))
for (i in 1:3)
  fit <- nlminb(model$env$last.par.best, model$fn, model$gr)

NLL <- fit$obj #objective value
pars <- length(fit$par)
-2*NLL + 2*pars + (2*pars*(pars+1))/(nyears-pars-1)

# Period = 4
map <- list(log_period = factor(rep(NA, njourn)))
model <- MakeADFun(mod_data, parameters, DLL="ex3", map = map,
                   silent = T, hessian = T)

fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
for (i in 1:3)
  fit <- nlminb(model$env$last.par.best, model$fn, model$gr)

NLL2 <- fit$obj #objective value
pars2 <- length(fit$par)
-2*NLL2 + 2*pars2 # period fixed at 4
-2*NLL + 2*pars # period estimated

# Period fixed at 4 is the best fit.

# ------
mod_data <- list(model = 2, nyears = nyears, njourn = njourn, 
                 data = as.matrix(data)[,1:2])

#parameter section
parameters <- list(b0 = rep(15, njourn), log_b1 = rep(0.02, njourn), 
                   amp = rep(0.2, njourn), phase = rep(0, njourn), 
                   log_period = rep(log(4), njourn))

# map <- list(log_period = factor(rep(NA, njourn)))

model <- MakeADFun(mod_data, parameters, DLL="ex3",
                   # map = map,
                   silent = T, hessian = T)

fit <- nlminb(model$par, model$fn, model$gr,
              control=list(eval.max=100000,iter.max=1000))
for (i in 1:3)
  fit <- nlminb(model$env$last.par.best, model$fn, model$gr)

# best <- model$env$last.par.best
# print(best)
# rep <- sdreport(model)
# print(rep)
# VarCo <- solve(model$he())
# # Check for Hessian
# print(sqrt(diag(VarCo)))

pred <- model$report()$cite_hat

plot(data$year, data$google, type = 'o')
lines(data$year, pred[,1], type = 'l',
      col = "blue")

plot(data$year, data$web, type = 'o')
lines(data$year, pred[,2], type = 'l',
      col = "blue")
