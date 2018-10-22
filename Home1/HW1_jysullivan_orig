setwd("Home1")

library(tidyverse); library(stats4); library(lme4)

# Question 1 ----

# 6 streams, 3 data points within each stream 
dat <- read.table("HOME1A.txt")

# Per Andre's email, only use 1st col (stream) and 3rd col (data)
dat <- dat[c(1,3)]

colnames(dat) <- c("stream","y")
N <- length(unique(dat$stream))
dat$replicate <- rep(1:3, N)

# fit model with lme4
fit <- lmer(formula = y ~ 1 + (1 | stream), data = dat)
summary(fit)

# beta = 78; sigma_b = 1; sigma_r = 2.5 - I tested a lot of different values and
# was having the issue of getting NaNs. 

# Problem 1
f <- function(beta, sigma_b, sigma_r) {
 
   require(dplyr)

  # sigma_r <- ifelse(sigma_r < 0, 1E-5, sigma_r) ## coerce to very small if neg
  # log_sigma_r <- log(sigma_r)
  # # log_sigma_b <- log(sigma_b)  
  
  nll <- 0 # initialize

  h <- 0.1 # stepsize to approximate integral in respect to b (random effects)
  bvec <- seq(-5, 5, h) # vector of b values
  simpson <-  h/3 * c(1, rep(c(4,2), length(bvec)/2 - 1), 4, 1) # lect2, slide 10 simpson's rule
  
  b <- dnorm(bvec, 0, 1) # random effects b_i~N(0,1), first part of Eqn 5
  
  tmp <- vector() # placeholder for each stream's integral
    
  for(i in 1:length(unique(dat$stream))) {  
    for(x in 1:length(bvec)) {  
      kprod <- 1 # placeholder
      for(k in 1:length(unique(dat$replicate))) {
        # Extract stream density obs
        y_ik <- filter(dat, stream==i & replicate==k) %>% pull(y)
        # Get likelihood for the data given mean, random effect, and variance
        # between streams (sigma_r) and within a stream (sigma_b), eqn 3 and
        # second half of eqn 5
        mu <- beta + b[x] * sigma_b 
        like <-  dnorm(y_ik, mu, sigma_r)
        kprod <- kprod * like  # product over k sites
      }
    
    tmp[x] <- simpson[x] * b[x] * kprod 
    nll <- nll - log(sum(tmp)) 
    }
  }
  return(nll)
}

# myfit <- stats4::mle(minuslogl = f, start = list(beta = 78, sigma_b = 16, sigma_r = 6),
#                      method = "L-BFGS-B", lower = c(60, 2.5, 3.5),
#                      upper = c(90, 20, 8))

myfit <- stats4::mle(minuslogl = f, start = list(beta = 78, sigma_b = 16, sigma_r = 6),
                     method = "BFGS")
# 
# Coefficients:
#   beta   sigma_b   sigma_r 
# 70.800368 57.542740  9.804852 

# Where beta is the mean, sigma_b is within stream variation, and sigma_r is the
# population or between stream bariation

# This took a really long time to run
# confint(myfit)
# 2.5 %   97.5 %
# beta    69.784281 71.79424 # did not contain the estimate from lme4
# sigma_b 53.990874 61.08924
# sigma_r  9.399112 10.23915

# GLARING OMISSION: No likelihood profile...

# Compare to lme4
summary(fit)
coef(fit)
ranef(fit)# beta 78.62 (std)

# what is the density of a typical
# stream? (ie mean population density) (fixed effect) 78.6

# what is the variance in density across streams (ie population variance)
# (fixed effect)? sigma_r^2 = 37.01

# What is the variation in density estimates within a stream? (random effect)
# 249.3

# append resids and fitted values
dat$resids <- resid(fit, type = "pearson")
dat$fitted <- fitted(fit)

# All residuals by fitted value
p <- ggplot(data = dat) +
  geom_point(aes(x = fitted, y = resids)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Pearsons residuals", x = "Fitted values") +
  theme_bw()
print(p) # looking good

# Residuals

# Residuals by stream - the assumption of homoskedasticity is met
p + facet_wrap(~stream)

# Extract random effects
randoms <- ranef(fit, condVar = TRUE)
qq <- attr(ranef(fit, condVar = TRUE)[[1]], "postVar")

df <- data.frame(Intercepts = randoms$stream[,1],
                 # error intervals (2 * sqrt(var))
                 sd.interc=2*sqrt(qq[,,1:length(qq)]),
                 lev.names=rownames(randoms$stream))

# Order output
df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])

p <- ggplot(data = df, 
            aes(lev.names, Intercepts)) +
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0, color="black") + 
  geom_point() +
  theme_bw() + 
  xlab("stream") + 
  ylab("") +
  coord_flip()

# Plot fixed effects

# plot data, X=fitted value
ggplot(data = dat, aes(y = y, x = factor(stream), colour = factor(stream))) +
  geom_point() +
  geom_point(data = dat, aes(y = fitted, x = factor(stream)), size = 4, shape = 4) +
  theme_bw() +
  labs(y = "Density", x = "Stream")

ggplot(data = dat, aes(x = log(ssb), colour = species)) +
  geom_point(aes(y = log(rec))) + 
  geom_line(aes(y = fitted, group = species)) +
  theme_bw()

# Errors are not not so normal.
qqnorm(dat2$resids, pch = 1, frame = FALSE)
qqline(dat2$resids, col = "steelblue", lwd = 2)

# Other attempt at Problem 1
# n <- 1000 # number of steps
# c <- 5 # upper bound
# a <- -5 # lower bound
# h <- (c-a)/n
# 
# b <- seq.int(a, c, length.out = n + 1)  
# out <- matrix(nrow = length(b), ncol = length(unique(dat$stream)))
# tmp <- vector(length = length(unique(dat$replicate)))
# 
# for(x in 1:length(b)) {
#   for(i in 1:length(unique(dat$stream))) {
#     for(k in 1:length(unique(dat$replicate))) {
#       y_ij <- filter(dat, stream==i & replicate==k) %>% pull(y)
#       # second part of eqn 5
#       tmp[k] <- 1/sqrt(2*pi*sigma_r) * exp((y_ij-beta-b[x]*sigma_b)^2/(2*sigma_r^2))
#     }
#     # Simpson's rule
#     if (b[x] == b[1] | x == b[n+1]) {
#       # eqn 5
#       out[x,i] <- (h/3) * (1/sqrt(2*pi) * exp(-(b[x]^2)/2) * prod(tmp))
#     } else {
#       if (b[x] %in% b[seq.int(2, length(b), 2)]) {
#         out[x,i] <- (h*4/3) * (1/sqrt(2*pi) * exp(-(b[x]^2)/2) * prod(tmp))
#       } else {
#         if (b[x] %in% b[seq.int(1, length(b), 2)]) {
#           out[x,i] <- (h*2/3) * (1/sqrt(2*pi) * exp(-(b[x]^2)/2) * prod(tmp))
#         }}}}}
# 
# colSums(out)
    
# Question 2 ----

dat2 <- read.csv("HOME1B.csv")
colnames(dat2) <- c("ssb", "rec", "species", "sbpr")
dat2$species <- as.factor(dat2$species)
str(dat2)

fit2 <- lmer(formula = log(rec) ~ log(ssb) - log(sbpr) - ssb + species - 1 + (1 | species), 
            data = dat2)

summary(fit2)

# append resids and fitted values
dat2$resids <- resid(fit2, type = "pearson")
dat2$fitted <- fitted(fit2)

# All residuals by fitted value
p <- ggplot(data = dat2) +
  geom_point(aes(x = fitted, y = resids)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(y = "Pearsons residuals", x = "Fitted values") +
  theme_bw()
print(p)# looks ok

# Residuals

# Residuals by species - residual patterns show heteroskedasticity (a violation
# of the constant variance assumption) some species (e.g. species 2, 3, 4, 7),
# indicating that variance structure may need some adjustment.
p + facet_wrap(~species)

# Extract random effects
randoms <- ranef(fit2, condVar = TRUE)
qq <- attr(ranef(fit2, condVar = TRUE)[[1]], "postVar")

df <- data.frame(Intercepts = randoms$species[,1],
           # error intervals (2 * sqrt(var))
           sd.interc=2*sqrt(qq[,,1:length(qq)]),
           lev.names=rownames(randoms$species))

# Order output
df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])

p <- ggplot(data = df, 
            aes(lev.names, Intercepts)) +
  geom_hline(yintercept=0) +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0, color="black") + 
  geom_point() +
  theme_bw() + 
  xlab("Species") + 
  ylab("") +
  coord_flip()

print(p)

# Plot fixed effects
ggplot(data = dat2, aes(x = log(ssb), colour = species)) +
  geom_point(aes(y = log(rec))) + 
  geom_line(aes(y = fitted, group = species)) +
  theme_bw()

# Errors are not not so normal.
qqnorm(dat2$resids, pch = 1, frame = FALSE)
qqline(dat2$resids, col = "steelblue", lwd = 2)
