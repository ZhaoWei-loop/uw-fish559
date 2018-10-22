# Homework 2
# Jane Sullivan
# jane.sullivan1@alaska.gov
# 2018-10-01

library(tidyverse)
library(tidyr)
library(TMB)
library(cowplot)
library(knitr)
source("tools.r")
setwd("Home2")

# Prep data 
Nyear <- readVec(file = "Home2.dat", 
                  string = "# Number of data points")
dat <- readMat(file = "Home2.dat", 
               string = "#Year M F", nrow = Nyear)

colnames(dat) <- c("Year", "M", "F")
dat <- as.data.frame(dat)

dat %>% 
  filter(M > 0 & F > 0) %>% 
  mutate(pM = M / (M + F) ) -> dat

with(dat, plot(x = Year, y = pM, type = 'o'))

data <- list(MM = dat$M, FF = dat$F)
parameters <- list(p = 0.4, tau = 1.05, dummy=2)

# 2. Fit model
# Compile and link
compile("HW2_jysullivan.cpp")
dyn.load(dynlib("HW2_jysullivan"))

# Model 1 ----
data$model <- 1
map <- list(dummy = factor(NA), tau = factor(NA))
model <- MakeADFun(data, parameters, map=map, 
                   DLL="HW2_jysullivan",silent=T,
                   hessian=T)
# checking for minimization
xx <- model$fn(model$env$last.par)
print(model$report())
fit <- nlminb(model$par, model$fn, model$gr, 
              control=list(eval.max=1000000,iter.max=100000))
best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(rep)
cat(model$report()$pfit,"\n")
res <- data.frame(p = model$report()$p,
           p_est = rep(as.list(rep, what = "Estimate")$`p`, length(model$report()$p)),
           p_std = rep(as.list(rep, what = "Std")$`p`, length(model$report()$p)))

ggplot(data = res, aes(x = p)) +
  geom_histogram(bins = 20, alpha = 0.5, aes(x = p, y = ..density..), position = "identity") +
  geom_vline(xintercept = res$p_est) +
  geom_vline(xintercept = res$p_est + (1.96 * res$p_std), lty = 3) +
  geom_vline(xintercept = res$p_est - (1.96 * res$p_std), lty = 3) +
  ggtitle(label = "Model 1") +
  theme(plot.title = element_text(hjust = 0)) -> mod1

ggsave("mod1_fit.png", plot = mod1, dpi = 300, height = 2.5, width = 5, units = "in")

# 4b. Model 2 ----

data$model <- 2
map <- list(dummy = factor(NA))
model <- MakeADFun(data, parameters, map=map, 
                   DLL="HW2_jysullivan",silent=T,
                   hessian=T)
# checking for minimization
xx <- model$fn(model$env$last.par)
print(model$report())
# Bounds on the parameters - needs to be equal to the number of ESTIMATED
# parameters not equal to the length of the total number of pars in the
# parameters vectors (dummy disappears) use:
length(model$par) # how many parameters for upper and lower bounds
lowbnd= c(0, # p
          0) # tau

uppbnd= c(1, # p
          Inf) # tau

fit <- nlminb(model$par, model$fn, model$gr, 
              control=list(rel.tol=1e-12,
                           eval.max=100000,iter.max=10000),
              lower=lowbnd,upper=uppbnd)

best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(rep)
cat(model$report()$pfit,"\n")
cat(model$report()$log_sigma_y,"\n")
model$report()$p
as.list(rep, what = "Std")$`p`

res <- data.frame(p = model$report()$p,
                  p_est = rep(as.list(rep, what = "Estimate")$`p`, length(model$report()$p)),
                  p_std = rep(as.list(rep, what = "Std")$`p`, length(model$report()$p)))

ggplot(data = res, aes(x = p)) +
  geom_histogram(bins = 20, alpha = 0.5, aes(x = p, y = ..density..), position = "identity") +
  geom_vline(xintercept = res$p_est) +
  geom_vline(xintercept = res$p_est + (1.96 * res$p_std), lty = 3) +
  geom_vline(xintercept = res$p_est - (1.96 * res$p_std), lty = 3) +
  ggtitle(label = "Model 2") +
  theme(plot.title = element_text(hjust = 0)) -> mod2

ggsave("mod2_fit.png", plot = mod2, dpi = 300, height = 2.5, width = 5, units = "in")

# 3. Simulation ----

set.seed(222)
sim_sexratio <- function(reps = 1000, # number of simulations
                         nyr = 25, # number of years 
                         nsamps = 100, # number of samples collected each yr)
                         mod = 1) { # which model to use
  
  # True p (expected value for p when p ~ beta(2,1))
  p_true <- mean(rbeta(reps*nyr, 2, 1))
  
  # Different p for each year and simulation 
  dat_ls <- list()
  MM <- vector()
  FF <- vector()
  
  for(i in 1:reps){
    # simulated p's, where p~beta(2,1)
    yr_p <- rbeta(nyr, 2, 1)
    
    for (j in 1:nyr) {
      # simulate data using the binomial distribution and our known p
      MM[j] <- rbinom(n = 1, size = nsamps, prob = yr_p[j])
      FF[j] <- 100 - MM[j]
    }
    dat_ls[[i]] <- list(MM = MM, FF = FF)
  }
  
  results <- list()
  bool <- vector()
  year <- vector()
  
  # mod = 2
  for(i in 1:reps) {
    
    data <- dat_ls[[i]]
    
    if(mod == 1) { 
      data$model <- 1
    } else {
      data$model <- 2
    } 
    
    model <- MakeADFun(data, parameters, map=map, 
                       DLL="HW2_jysullivan",silent=T,hessian=T)
    
    fit <- nlminb(model$par, model$fn, model$gr, 
                  control=list(rel.tol=1e-12,
                               eval.max=100000,iter.max=10000),
                  lower=lowbnd,upper=uppbnd)
    
    best <- model$env$last.par.best
    rep <- sdreport(model)
    
    p_est <- as.list(rep, what = "Est")$`p`
    p_se <- as.list(rep, what = "Std")$`p`
    p_lower95 <- p_est - 1.96 * p_se
    p_upper95 <- p_est + 1.96 * p_se
    
    for(j in 1:nyr) {
      bool[j] <- ifelse(p_true <= p_upper95 & p_true >= p_lower95, 1, 0)
      year[j] <- j
    }
    
    results[[i]] <- data.frame(p_est = rep(p_est, nyr),
                               p_se = rep(p_se, nyr),
                               p_lower95 = rep(p_lower95, nyr),
                               p_upper95 = rep(p_upper95, nyr),
                               p_true = rep(p_true, nyr), bool = bool,
                               year = year, sim = i)
    
  }
  return(results)
}

# Sim Model 1 ----
lowbnd= c(0) # p
uppbnd= c(1) # p
map <- list(dummy = factor(NA), tau = factor(NA))
results <- sim_sexratio(reps = 1000, # number of simulations
                        nyr = 25, # number of years 
                        nsamps = 100, # number of samples collected each yr)
                        mod = 1)
res <- do.call("rbind", results)

res %>% 
  group_by(bool) %>% 
  summarize(Count = n()) %>% 
  mutate(bool = ifelse(bool==0, "No", "Yes"),
         Percent = formatC(round(Count/length(res$p_est)*100, 1), small.interval = 1)) %>% 
  select(`95% CI contains true p?` = bool, Count, Percent) %>% 
  kable()

res %>% 
  ggplot(aes(x = p_est, fill = factor(bool), colour = factor(bool))) + 
  geom_histogram(bins = 100, alpha = 0.5, 
                 aes(y = ..density.., fill = factor(bool), colour = factor(bool)), 
                 position = "identity") +
  geom_vline(xintercept = res$p_true, size = 1, lty = 2) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  scale_colour_grey(start = 0.3, end = 0.8) +
  theme(legend.position = "none") + 
  ggtitle(label = "Model 1 simulation") +
  labs(x = "Estimated p", y = "Density") +
  theme(plot.title = element_text(hjust = 0)) -> psims

ggsave("psims_model1.png", plot = psims, dpi = 300, height = 2.5, width = 5, units = "in")

# Sim Model 2 ----
map <- list(dummy = factor(NA))
lowbnd= c(0, # p
          0) # tau

uppbnd= c(1, # p
          Inf) # tau
results <- sim_sexratio(reps = 1000, # number of simulations
                        nyr = 25, # number of years 
                        nsamps = 100, # number of samples collected each yr)
                        mod = 2)
res <- do.call("rbind", results)

res %>% 
  group_by(bool) %>% 
  summarize(Count = n()) %>% 
  mutate(bool = ifelse(bool==0, "No", "Yes"),
         Percent = formatC(round(Count/length(res$p_est)*100, 1), small.interval = 1)) %>% 
  select(`95% CI contains true p?` = bool, Count, Percent) %>% 
  kable()

res %>% 
  ggplot(aes(x = p_est, fill = factor(bool), colour = factor(bool))) + 
  geom_histogram(bins = 70, alpha = 0.5, 
                 aes(y = ..density.., fill = factor(bool), colour = factor(bool)), 
                 position = "identity") +
  geom_vline(xintercept = res$p_true, size = 1, lty = 2) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  scale_colour_grey(start = 0.3, end = 0.8) +
  theme(legend.position = "none") + 
  ggtitle(label = "Model 2 simulation") +
  labs(x = "Estimated p", y = "Density") +
  theme(plot.title = element_text(hjust = 0))-> psims

ggsave("psims_model2.png", plot = psims, dpi = 300, height = 2.5, width = 5, units = "in")

dat
ggplot(dat, aes(x = Year, y = pM)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  ylab("Proportion male") -> ts

ggsave("ts.png", plot = ts, dpi = 300, height = 2.5, width = 5, units = "in")

acf(dat$pM)
