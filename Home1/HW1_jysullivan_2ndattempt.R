setwd("Home1")

library(tidyverse); library(stats4); library(lme4); library(nlme)

# 1a ----

# 6 streams, 3 data points within each stream 
dat <- read.table("HOME1A.txt")

# Per Andre's email, only use 1st col (stream) and 3rd col (data). Each row is a
# stream, each column is a replicate
dat <- matrix(dat[,3], ncol = 3, byrow = TRUE)

h <- 0.025 # stepsize to approximate integral in respect to b (random effects)
bvec <- seq(-10, 10, h) # vector of b values
simpson <-  h/3 * c(1, rep(c(4,2), length(bvec)/2 - 1), 4, 1) # lect2, slide 10 simpson's rule
  
f <- function(beta, sigma_b, sigma_r, dat) {
 
   #require(dplyr)

  nll <- 0 # initialize

  b <- dnorm(bvec, 0, 1) # random effects b_i~N(0,1), first part of Eqn 5
  
  # tmp <- vector() # placeholder for each stream's integral
    
  for(i in 1:length(dat[,1])) { # loop over each stream
    
    kprod <- rep(1, length(b)) # placeholder
    
    for(k in 1:length(dat[1,])) { # loop over each replicate in a stream
      
      # Get likelihood for the data given mean, random effect, and variance
      # between streams (sigma_r) and within a stream (sigma_b), eqn 3 and
      # second half of eqn 5
      like <- dnorm(dat[i,k], (beta + bvec * sigma_b), sigma_r)
      kprod <- kprod * like  # product over k sites
    }
    
    kprod <- kprod * b 
    tmp <- sum(simpson * kprod)
    nll <- nll - log(tmp) # negative log likelihood
    
  }
  return(nll)
}

# Test function
fixed_pars <- c(70, 10, 6)
tst <- f(fixed_pars[1], fixed_pars[2], fixed_pars[3], dat)
print(tst)

myfit <- stats4::mle(minuslogl = f, start = list(beta = 70, sigma_b = 10, sigma_r = 6),
                   fixed = list(dat = dat), method = "BFGS")

summary(myfit)
(mle <- -1*logLik(myfit)[1]) # maximum likelihood
coef(myfit) # parameter estimates and fitted vals

# Where beta is the mean, sigma_b is within stream variation, and sigma_r is the
# population or between stream bariation

# 1b ----

# Likelihood profile

# First run it with the simple/coarse vector of values
beta_vec <- seq(55, 105, by = 1)
# Then add high resolution to the areas close to the interval
beta_vec <- c(beta_vec, seq(from=63,to=65,by=0.01),seq(from=91,to=93,by=0.01))
beta_vec <- sort(beta_vec)

like_prof <- NULL

for (i in beta_vec) {
  
  # Fix beta to one of the beta values in beta_vec. Note using the Nelder-Mead
  # algorithm
  f1 <- mle(f, start = list(sigma_b = 16, sigma_r = 6),
            fixed = list(beta = i, dat = dat), 
            method = "Nelder-Mead")
  
  # Add to the vector of negative log likelihoods
  like_prof <- c(like_prof, -1*logLik(f1)[1])
}

data.frame(Beta = beta_vec, NLL = like_prof - mle) -> df

# We know that the lower and upper bounds are somewhere between 64 and 65 and 92
# and 93 respectively
df %>% 
  mutate(nll = round(NLL, 2)) %>% 
  filter(nll == 1.92) 

# Use uniroot to find the exact value
ci_intervals <- list(lower = c(64, 65),
                     upper = c(92, 93))

# Adapted JBest's code to mine to get CI using uniroot(). He used 2*NLL and
# 2*192, but it worked fine without that
ci_uniroots <- lapply(ci_intervals,
                      function(interval) {
                        uniroot(function(beta) {

                          f1 <- mle(f, start = list(sigma_b = 16, sigma_r = 6),
                              fixed = list(beta = beta, dat = dat), 
                              method = "Nelder-Mead")
                          
                          (-1*logLik(f1)[1] - mle) - 1.92 # qchisq(0.95, 1) / 2 = 1.92
                          
                        }, interval = interval)})
ci <- lapply(ci_uniroots,
             function(root) root$root)

print(ci)

df %>% 
  ggplot(aes(x = Beta, y = NLL)) +
  geom_line() +
  geom_hline(yintercept = 1.92, lty = 2, col = "grey") +
  geom_vline(xintercept = ci$lower[1], lty = 2, col = "grey") + 
  geom_vline(xintercept = ci$upper[1], lty = 2, col = "grey") -> lprof

ggsave("likeprof.png", plot = lprof, dpi = 300, height = 2.5, width = 5, units = "in")

# 1c - lme4 ----

# fit model with lme4
dat <- read.table("HOME1A.txt")
dat <- dat[, c(1,3)] # don't need 2nd column
colnames(dat) <- c("stream", "y")
dat$stream <- as.factor(dat$stream)

fit_lme4 <- lmer(formula = y ~ 1 + (1 | stream), data = dat)
summary(fit_lme4)
confint(fit_lme4)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# stream   (Intercept) 249.29   15.789  
# Residual              37.01    6.084  
# Number of obs: 18, groups:  stream, 6
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   78.629      6.603   11.91

# what is the density of a typical
# stream? (ie mean population density) (fixed effect) 78.6

# what is the variance in density across streams (ie population variance)
# (fixed effect)? sigma_r^2 = 37.01

# What is the variation in density estimates within a stream? (random effect)
# 249.3

# append resids and fitted values
dat$resids <- resid(fit_lme4, type = "pearson")
dat$fitted <- fitted(fit_lme4)

# All residuals by fitted value
p <- ggplot(data = dat) +
  geom_point(aes(x = fitted, y = resids)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Pearsons residuals", x = "Fitted values") +
  theme_bw()
print(p) # looking good

ggsave("resid_allstream.png", plot = p, dpi = 300, height = 2.5, width = 5, units = "in")

# Residuals by stream - the assumption of homoskedasticity is met
p + facet_wrap(~stream) -> p
print(p)

ggsave("resid_bystream.png", plot = p, dpi = 300, height = 5, width = 5, units = "in")

# Extract random effects
randoms <- ranef(fit_lme4, condVar = TRUE)
qq <- attr(ranef(fit_lme4, condVar = TRUE)[[1]], "postVar")

df <- data.frame(Intercepts = randoms$stream[,1],
                 # error intervals (2 * sqrt(var))
                 sd.interc=2*sqrt(qq[,,1:length(qq)]),
                 lev.names=rownames(randoms$stream))

# Order output
df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])

p <- ggplot(data = df, 
            aes(lev.names, Intercepts)) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey") +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0, color="black") + 
  geom_point() +
  theme_bw() + 
  xlab("stream") + 
  ylab("") +
  coord_flip()

print(p)

ggsave("stream_re.png", plot = p, dpi = 300, height = 2.5, width = 5, units = "in")

# Plot fixed effects

# plot data, X=fitted value
ggplot(data = dat, aes(y = y, x = stream)) +
  geom_point() +
  geom_point(data = dat, aes(y = fitted, x = factor(stream)), size = 4, shape = 4) +
  theme_bw() +
  labs(y = "Density", x = "Stream") -> p
print(p)
ggsave("streams_fixedeff.png", plot = p, dpi = 300, height = 2.5, width = 5, units = "in")

# Errors are so so normal.
qqnorm(dat$resids, pch = 1, frame = FALSE)
qqline(dat$resids, col = "black", lwd = 2)

# 1c - nlme
# For comparison, fit the same model with nlme::lme

# necessary step for nlme
dat_nlme <- groupedData(y ~ 1 | stream, data = dat)
fit_nlme <- lme(fixed = y ~ 1, random = ~ 1 | stream, data = dat, method = "ML") 
summary(fit_nlme)
coef(fit_nlme)
# 2 lme4

dat2 <- read.csv("HOME1B.csv")
colnames(dat2) <- c("ssb", "rec", "species", "sbpr")
dat2$species <- as.factor(dat2$species)
dat2$y <- log(dat2$rec * dat2$sbpr / dat2$ssb)

str(dat2)

fit2_lme4 <- lmer(formula = y ~ species : ssb - 1 + (1 | species), 
            data = dat2)

summary(fit2_lme4)

# append resids and fitted values
dat2$resids <- resid(fit2_lme4, type = "pearson")
dat2$fitted <- fitted(fit2_lme4)

# All residuals by fitted value
p <- ggplot(data = dat2) +
  geom_point(aes(x = fitted, y = resids)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Pearsons residuals", x = "Fitted values") +
  theme_bw()
print(p) # looks like there is some heteroskedasticity

ggsave("resids_sr.png", plot = p, dpi = 300, height = 2.5, width = 5, units = "in")

# Residuals by species - residual patterns show heteroskedasticity (a violation
# of the constant variance assumption) some species (e.g. species 2, 3, 4, 7),
# indicating that variance structure may need some adjustment.
p + facet_wrap(~species) -> p
print(p)

ggsave("resids_byspp.png", plot = p, dpi = 300, height = 5, width = 5, units = "in")

# Extract random effects
randoms <- ranef(fit2_lme4, condVar = TRUE)
qq <- attr(ranef(fit2_lme4, condVar = TRUE)[[1]], "postVar")

df <- data.frame(Intercepts = randoms$species[,1],
           # error intervals (2 * sqrt(var))
           sd.interc=2*sqrt(qq[,,1:length(qq)]),
           lev.names=rownames(randoms$species))

# Order output
df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])

p <- ggplot(data = df, 
            aes(lev.names, Intercepts)) +
  geom_hline(yintercept=0, lty = 2, colour = "grey") +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0, color="black") +
  geom_point() +
  theme_bw() + 
  xlab("Species") + 
  ylab("") +
  coord_flip()

print(p)

ggsave("species_re.png", plot = p, dpi = 300, height = 2.5, width = 5, units = "in")


# Plot fixed effects
p <- ggplot(data = dat2, aes(x = ssb)) +
  geom_point(aes(y = y)) + 
  geom_line(aes(y = fitted, group = species), colour = "grey") +
  theme_bw() + 
  facet_wrap(~ species, scales = "free") 
print(p)

ggsave("spp_fixed_fitted.png", plot = p, dpi = 300, height = 5, width = 5, units = "in")

# Errors are not not so normal.
qqnorm(dat2$resids, pch = 1, frame = FALSE)
qqline(dat2$resids, col = "black", lwd = 2)
